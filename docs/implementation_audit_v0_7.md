# RegCompassR v0.6-v0.7 implementation audit

This audit checks the current tool code against `RegCompass_Multiome_adjusted_development_spec.md` for mathematical, statistical, and biological-computing consistency.

## Selected reaction demand QP planning

- `rc_select_reactions()` now implements the required selected set as the union of exchange reactions, transport reactions, top Layer 1 variable reactions, optional top sample-level differential reactions, and user-specified reactions.
- The variability component uses row-wise variance of the reaction-by-pool Layer 1 matrix `C_rel`; this is a screening statistic only and not a biological-replicate hypothesis test.
- Optional differential reaction input is deliberately external to `rc_select_reactions()` so it can come from sample-level statistics such as `rc_lm_by_reaction()` rather than direct pool-level tests.
- `rc_estimate_selected_demand_qp()` implements the required workload formula `N_QP = n_pools * (1 + n_selected_reactions)` and reports serial/parallel time estimates plus checkpoint count.

## Sample-level statistics

- `rc_sample_aggregate()` uses median aggregation over pools within each sample × annotated cell type group. This matches the planned default statistical unit and avoids treating pools as independent biological replicates.
- `rc_lm_by_reaction()` fits ordinary sample-level linear models reaction by reaction. This matches the v0.7 plan to start with simple sample-level linear models and defer mixed models.
- P-values are BH-adjusted within each model term across reactions, consistent with reaction-wise multiple testing within a contrast.

## Regulator ranking

- `rc_rank_regulators()` is candidate prioritization only, not causal driver inference.
- Ranking is performed within each reaction, so regulator candidates for unrelated reactions are not compared on a single global scale.
- Multiple evidence columns are converted to ranks and combined using an order-statistic robust rank aggregation score. A BH q-value is then computed within each reaction across candidate regulators.
- Motif and enhancer columns are interpreted as support tiers rather than proof of reaction activation or causality.

## Execution constraints

- Long selected-demand QP runs expose `BPPARAM` for Linux parallelism.
- Serial checkpoint/resume is supported for selected-demand QP sweeps. Shared checkpoint files are intentionally disallowed with `BPPARAM` to avoid unsafe concurrent writes.
