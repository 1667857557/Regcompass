# Pando and multitask GRN parameter policy

This document is the canonical parameter contract for `grn_mode = "multitask_shared_backbone"` in RegCompassR 1.8.10.

## Design invariants

1. Every condition of one cell type uses the same structural Pando candidate dictionary.
2. Candidate removal before fitting cannot use full-data target correlation, fitted effect size, or P values.
3. Condition coefficients are estimated on one shared edge scale.
4. Elastic-net sparsity acts directly on the condition coefficient \(\theta_{e,c}\), permitting exact condition-specific zeros.
5. The reported global backbone is a derived cross-condition mean, not a separately penalised latent coefficient.
6. Active targets must improve on a condition-centred out-of-fold null model.
7. Bootstrap resampling uses biological sample/donor clusters when a valid `sample_col` is available and records any fallback explicitly.

## Canonical defaults

### Pando structural design

| Parameter | Default | Role |
|---|---:|---|
| `min_cells` | `100` | Minimum cells in every condition of a cell type. |
| `peak_to_gene_method` | `"GREAT"` | GREAT-style basal-plus-extension regulatory-domain rule. |
| `upstream` | `100000` | Basal regulatory domain upstream of the gene/TSS anchor. |
| `downstream` | `0` | Basal regulatory domain downstream of the gene/TSS anchor. |
| `extend` | `1000000` | Maximum GREAT domain extension. |
| `only_tss` | `FALSE` | Use the annotated gene body rather than TSS-only anchoring. |
| `min_tf_detection` | `0` | Preserve condition-restricted TFs in pooled preprocessing. |
| `min_peak_detection` | `0` | Preserve condition-restricted peaks in pooled preprocessing. |
| `min_target_detection` | `0` | Preserve condition-restricted targets in pooled preprocessing. |
| `max_edges_per_target` | `Inf` | Preserve the complete structural candidate set. |

`GREAT` is canonical because it allows distal structural hypotheses while keeping candidate construction independent of target-expression correlation and fitted effects. The condition-aware observability filter and direct-theta regularisation must remain enabled. `Signac` remains an explicit sensitivity option.

### Multitask model

| Parameter | Default | Role |
|---|---:|---|
| `alpha` | `0.5` | L1 sparsity plus ridge stabilisation. |
| `global_penalty_factor` | `1` | Compatibility alias for the common direct-theta penalty. |
| `deviation_penalty_factor` | `1` | Compatibility alias; must equal `global_penalty_factor`. |
| `lambda_rule` | `"lambda.1se"` | Conservative lambda from leakage-resistant CV. |
| `nfolds` | `5` | Condition-stratified cell-level folds. |
| `n_bootstrap` | `100` | Condition-stratified bootstrap replicates. |
| `min_selection_frequency` | `0.7` | Minimum direct-theta bootstrap selection frequency. |
| `min_sign_stability` | `0.8` | Minimum conditional sign stability. |
| `min_bootstrap_success_fraction` | `0.8` | Minimum completed bootstrap fraction. |
| `min_cv_rsq` | `0` | User floor; activation still requires positive OOF R-squared. |
| `candidate_screen_threshold` | `0` | Full-data outcome screening is disabled. |
| `max_edges_per_target` | `Inf` | No deterministic top-K truncation. |
| `min_detected_cells_per_condition` | `10` | Absolute observability floor. |
| `min_detection_fraction_per_condition` | `0.01` | Relative observability floor. |

### Bootstrap sample metadata

| Argument | Default | Role |
|---|---:|---|
| `sample_col` | `NULL` | Biological sample/donor column for Stage 1 cluster bootstrap. |

A valid `sample_col` activates sample-cluster bootstrap. `sample_col = NULL` or a named column that does not exist prints a warning with the exact reason and falls back to condition-stratified cell resampling. An existing sample column with missing or empty IDs is rejected. Fewer than two samples in a condition produces a low-replication warning but does not silently switch resampling units.

`sample_col` is not a Stage 2 metacell stratum and does not create sample-specific metabolic models.

## Shared structural and observable universes

Pando returns

\[
\mathcal U_m^{struct}
=\{e=(t,p,g):\text{motif and GREAT domains support }e\}
\]

for cell type \(m\), using all conditions together.

For condition \(c\), define

\[
m_c=\min\left(n_c,\max\left(10,\lceil0.01n_c\rceil\right)\right).
\]

The shared model universe is

\[
\mathcal U_m^{model}
=\left\{e\in\mathcal U_m^{struct}:\exists c,
N^{TF\times peak}_{e,c}\ge m_c,
N^{target}_{e,c}\ge m_c\right\}.
\]

The rule uses detection only. Structural and model universes receive separate MD5 fingerprints.

## Direct condition-theta objective

For predictor

\[
x_{e,u}=T_{t,u}A_{p,u},
\]

target and predictors are centred within condition and predictors are divided by one edge scale shared across conditions:

\[
y_u^\circ=\sum_e\widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

The weighted elastic-net objective is

\[
\begin{aligned}
\min_\Theta\quad
&\frac12\sum_u w_u
\left(y_u^\circ-\sum_e\widetilde x_{e,u}\theta_{e,c(u)}\right)^2\\
&+\lambda\alpha p_\theta\sum_{e,c}|\theta_{e,c}|\\
&+\frac{\lambda(1-\alpha)p_\theta}{2}\sum_{e,c}\theta_{e,c}^2,
\end{aligned}
\]

where \(w_u\propto1/n_{c(u)}\).

The reported summaries are

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

`global_penalty_factor` and `deviation_penalty_factor` refer to the same \(p_\theta\) and must be equal.

## Cross-validation

Five folds are assigned separately within every condition. For each fold:

1. target and predictor condition means are estimated on training cells;
2. training means are applied to validation cells;
3. predictor scales are estimated on training cells;
4. training scales are applied to validation cells;
5. every lambda is evaluated on validation residuals.

The selected lambda uses the one-standard-error rule. Active edges require strictly positive joint target-model out-of-fold \(R^2\). Cross-validation remains cell-level and is not a biological-replicate significance test.

## Bootstrap reproducibility

For condition \(c\), let \(D_c\) denote its observed sample IDs.

### Sample-cluster mode

A replicate draws \(|D_c|\) sample IDs with replacement from \(D_c\). Every selected sample contributes all cells belonging to that sample and condition. Donor cluster sizes are preserved, so replicate cell counts may vary.

### Cell fallback mode

A replicate draws \(n_c\) cells with replacement from condition \(c\). This is used only when `sample_col` is not supplied or the named column does not exist, and the fallback reason is printed and recorded.

Both modes re-centre the resampled matrices and refit at the selected lambda. Selection frequency is

\[
\widehat\Pi_{e,c}
=\frac1{B_s}\sum_{b\in\mathcal B_s}
I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon).
\]

Conditional sign stability is

\[
\rho_{e,c}=\left|
\frac{\sum_bI^{(b)}_{e,c}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_bI^{(b)}_{e,c}}
\right|.
\]

The reliability-weighted estimate is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\widehat\Pi_{e,c}\rho_{e,c}.
\]

Sample-cluster bootstrap measures reproducibility across observed sample clusters, not a treatment effect. Cell fallback measures cell-resampling reproducibility only. Neither substitutes for donor-level differential inference.

## Active-edge rule

An edge is active only when all conditions hold:

\[
\begin{aligned}
&B_s/B\ge0.8,\\
&\widehat\Pi_{e,c}\ge0.7,\\
&\rho_{e,c}\ge0.8,\\
&|\widehat\theta_{e,c}|>\max(\tau_{effect},\varepsilon),\\
&R^2_{CV}>0,\\
&R^2_{CV}\ge\tau_{CV}\quad\text{when }\tau_{CV}>0.
\end{aligned}
\]

`padj` remains `NA` because this is bootstrap stability selection rather than a classical per-edge hypothesis test.

## Auditable output fields

Stage 1 reports:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
selection_frequency
selection_frequency_mc_se
selection_frequency_lower_95
selection_frequency_upper_95
sign_stability
```

The compact final result carries the selected bootstrap policy in result parameters and `stage_provenance`. Timing is printed in R and is not persisted in result objects or files.
