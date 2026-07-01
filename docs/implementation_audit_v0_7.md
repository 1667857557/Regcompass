# RegCompassR v0.1-v0.7 implementation audit

This audit checks the current tool code against `RegCompass_Multiome_adjusted_development_spec.md` for mathematical, statistical, and biological-computing consistency. The package is still an early milestone implementation: modules marked **implemented MVP** are expected to run, while modules marked **planned / not implemented** must not be interpreted as available analysis functionality.

## Module-by-module status

| Spec module | Current code status | Mathematical / biological audit |
|---|---|---|
| Module 0: QC, normalization, joint embedding | **Input validation only** via Seurat object checks. | Correctly assumes user input is an already annotated Seurat/Signac multiome object and does not rerun clustering, WNN, normalization, or embedding. |
| Module 1: sample-aware micropooling | **Implemented MVP**. | Pools are formed within sample-aware strata so cells from different biological samples are not mixed. This preserves valid downstream sample-level statistics. |
| Module 2: gene-level multiome support | **Not fully implemented**. | Current code does not claim complete enhancer/motif/peak-gene support modeling; ATAC should remain support/confidence rather than a flux bound. |
| Module 3: GPR-aware reaction capacity | **Implemented MVP** for simple GPR parsing, robust gene scores, Boltzmann AND, OR alternatives, Q95 calibration, and diagnostics. | Layer 1 capacity remains a relative reaction-capacity potential, not an inferred flux. |
| Module 4: GEM preprocessing / integrity | **Toy/minimal GEM validation only**. | Full Human-GEM preprocessing, flux consistency, reversibility splitting, compartment checks, and loop diagnostics remain outside current scope. |
| Module 5: network-constrained QP feasibility | **Toy/minimal QP MVP**. | Baseline QP uses sparse OSQP matrices and soft quadratic penalties, but is not a full Human-GEM production QP implementation. |
| Module 6: selected reaction-demand exact feasibility | **Implemented planning + toy selected-demand execution**. | Selection now covers top Layer 1 variable, optional top sample-level differential, exchange, transport, and user reactions. Workload estimation follows `N_QP = n_pools * (1 + n_selected)`. |
| Module 7: limited FVA | **Planned / not implemented**. | No production FVA should be claimed. |
| Module 8: pathway / meta-reaction aggregation | **Planned / not implemented**. | No pathway/meta-reaction scoring should be claimed yet. |
| Module 9: sample-aware statistics | **Implemented MVP**. | Pool-level scores are aggregated to sample × cell type medians before linear modeling, avoiding pseudo-replication from direct pool-level tests. |
| Module 10: signaling-anchored rewiring | **Planned / not implemented**. | No signaling-anchored causal or rewiring model should be claimed. |
| Module 11: regulator candidate prioritization | **Implemented MVP**. | Regulator rankings are candidate prioritization only. Evidence ranks are aggregated within each reaction and do not imply causal driver discovery. |
| Module 12: cell-wise projection | **Planned / not implemented**. | No cell-wise projection/interpolation should be claimed. |

## Selected reaction demand QP planning

- `rc_select_reactions()` implements the required selected set as the union of exchange reactions, transport reactions, top Layer 1 variable reactions, optional top sample-level differential reactions, and user-specified reactions.
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

## Important limitations to preserve scientific interpretation

1. Layer 1 reaction capacity is a relative potential score, not flux.
2. Selected demand-QP feasibility is exact only for selected reactions and the supplied model; it is not full reaction coverage.
3. ATAC/motif support is regulatory confidence/support and must not be used as proof of reaction activation.
4. Regulator ranking is candidate prioritization and must not be described as causal driver discovery.
5. Current QP/GEM support remains a minimal/toy MVP until full model preprocessing and integrity diagnostics are implemented.
