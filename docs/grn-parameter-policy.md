# Pando and multitask GRN parameter policy

This document is the canonical parameter contract for `grn_mode = "multitask_shared_backbone"` in RegCompassR 1.8.10.

## Design invariants

1. Every condition of one cell type uses the same structural Pando candidate dictionary.
2. Candidate removal before fitting cannot use full-data target correlation, fitted effect size or P values.
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
| `min_tf_detection` | `0` | Do not remove condition-restricted TFs in pooled Pando preprocessing. |
| `min_peak_detection` | `0` | Do not remove condition-restricted peaks in pooled Pando preprocessing. |
| `min_target_detection` | `0` | Do not remove condition-restricted targets in pooled Pando preprocessing. |
| `max_edges_per_target` | `Inf` | Preserve the complete structural candidate set. |

`GREAT` is the canonical structural rule because it allows distal regulatory coverage while keeping candidate construction independent of target-expression correlation, fitted coefficients and condition effects. With `extend = 1000000`, it deliberately creates a broad hypothesis space; the condition-aware observability filter and direct-theta regularisation must therefore remain enabled. `Signac` remains an explicit sensitivity option, but is no longer the default.

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
| `min_cv_rsq` | `0` | User floor; activation still requires strictly positive out-of-fold R-squared. |
| `candidate_screen_threshold` | `0` | Full-data outcome screening is disabled. |
| `max_edges_per_target` | `Inf` | No deterministic top-K truncation. |
| `min_detected_cells_per_condition` | `10` | Absolute observability floor. |
| `min_detection_fraction_per_condition` | `0.01` | Relative observability floor. |

### Bootstrap sample metadata

| Argument | Default | Role |
|---|---:|---|
| `sample_col` | `NULL` | Biological sample/donor column for Stage 1 condition-stratified cluster bootstrap. |

A valid `sample_col` activates sample-cluster bootstrap. `sample_col = NULL` or a named column that does not exist prints a warning with the exact reason and falls back to condition-stratified cell resampling. An existing sample column with missing or empty IDs is rejected. Fewer than two samples in a condition produces a low-replication warning but does not silently switch resampling units.

`sample_col` is not a Stage 2 metacell stratum and does not create sample-specific metabolic models.

## Shared structural and observable universes

Pando returns

\[
\mathcal U_m^{struct}
=\{e=(t,p,g):\text{motif and GREAT peak-to-gene domains support }e\}
\]

for cell type \(m\), using all conditions together.

For condition \(c\), define

\[
m_c=\min\left(n_c,\max\left(10,\lceil0.01n_c\rceil\right)\right).
\]

For edge \(e=(t,p,g)\), let

\[
N^{TF\times peak}_{e,c}
=\sum_{u\in c}I(T_{t,u}>0\land A_{p,u}>0)
\]

and

\[
N^{target}_{e,c}=\sum_{u\in c}I(Y_{g,u}>0).
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

target and predictors are centred within condition and predictors are divided by one edge scale shared across conditions. The design contains one coefficient for every edge-condition pair:

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
&+\frac{\lambda(1-\alpha)p_\theta}{2}
\sum_{e,c}\theta_{e,c}^2,
\end{aligned}
\]

where

\[
w_u\propto\frac1{n_{c(u)}}.
\]

Because the L1 penalty is applied directly to \(\theta_{e,c}\), the solution can naturally contain

\[
\theta_{e,A}\ne0,\qquad\theta_{e,B}=0.
\]

The reported backbone and deviations are derived:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

This replaces the former parameterisation that penalised \((\beta,\gamma)\) and only indirectly produced sparse condition topology.

## Penalty aliases

For compatibility, `global_penalty_factor` and `deviation_penalty_factor` remain accepted. They refer to the same \(p_\theta\) and therefore must be equal. A configuration such as

```r
global_penalty_factor = 1
deviation_penalty_factor = 2
```

is rejected because the direct-theta model has no separate deviation-coordinate penalty. A stronger conserved-backbone prior would require a distinct fused or grouped multitask penalty and is not silently approximated.

## Cross-validation

Five folds are assigned separately within every condition. For each fold:

1. condition means for target and predictors are estimated on training cells;
2. training means are applied to validation cells;
3. predictor scales are estimated on training cells;
4. training scales are applied to validation cells;
5. every lambda is evaluated on validation residuals.

The selected lambda uses the one-standard-error rule. Active edges require strictly positive joint target-model out-of-fold \(R^2\). This is an internal predictive reliability statistic, not a biological-replicate significance test. Cross-validation remains cell-level even when bootstrap uses donor/sample clusters.

## Bootstrap reproducibility

For condition \(c\), let \(D_c\) denote its observed sample/donor IDs.

### Sample-cluster mode

A bootstrap replicate draws \(|D_c|\) sample IDs with replacement from \(D_c\). Every occurrence of a selected sample contributes all cells belonging to that sample and condition. Donor cluster sizes are preserved, so replicate cell counts may vary.

### Cell fallback mode

When `sample_col` is omitted or the named metadata column does not exist, RegCompass prints the exact fallback reason and draws \(n_c\) cells with replacement from condition \(c\). An existing column with missing or empty sample IDs is an error.

Both modes re-centre the resampled matrices and refit at the selected lambda. Selection frequency is based on direct condition coefficients:

\[
\widehat\Pi_{e,c}
=\frac1{B_s}\sum_{b\in\mathcal B_s}
I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon).
\]

Conditional sign stability is

\[
\rho_{e,c}=
\left|
\frac{\sum_bI^{(b)}_{e,c}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_bI^{(b)}_{e,c}}
\right|.
\]

The reliability-weighted estimate is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\widehat\Pi_{e,c}\rho_{e,c}.
\]

Sample-cluster bootstrap estimates reproducibility across the observed biological sample clusters. It does not provide a donor-level treatment-effect estimate or formal half-sample stability-selection error control. Cell fallback estimates cell-resampling reproducibility only.

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

## Output and timing contract

Stage 1 coefficient, stability, cell-type status, and condition-by-cell-type status outputs retain:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
```

The compact final result retains the workflow-level bootstrap policy and edge-level provenance before Stage 3 core construction. Each public stage prints elapsed time and final status in R; timing is not stored in stage objects, `result$timing`, `step_timing.tsv`, or `00_execution_timing.tsv`.
