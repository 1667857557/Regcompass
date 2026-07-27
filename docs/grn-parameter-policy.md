# Pando and multitask GRN parameter policy

This document is the canonical parameter contract for
`grn_mode = "multitask_shared_backbone"`.

## Design goals

The parameter system must preserve four invariants:

1. every condition of one cell type uses the same structural Pando candidate
   dictionary;
2. candidate removal before fitting cannot use full-data target correlation or
   fitted effect size;
3. condition-specific coefficients are estimated on a shared scale and must
   improve on an out-of-fold null predictor;
4. thresholds must have an interpretable statistical or biological role rather
   than encode an arbitrary top-K or fixed undocumented prior.

## Canonical defaults

### Pando structural design

| Parameter | Default | Role |
|---|---:|---|
| `min_cells` | `100` | Minimum cells in every condition of a cell type. |
| `peak_to_gene_method` | `"Signac"` | Canonical peak-to-gene rule supported by current Pando. |
| `upstream` | `100000` | Candidate region upstream of the gene. |
| `downstream` | `0` | Candidate region downstream of the gene. |
| `extend` | `1000000` | GREAT extension, used only when GREAT is selected. |
| `only_tss` | `FALSE` | Measure the Signac window from the gene body rather than TSS only. |
| `min_tf_detection` | `0` | Do not remove condition-restricted TFs in pooled Pando preprocessing. |
| `min_peak_detection` | `0` | Do not remove condition-restricted peaks in pooled Pando preprocessing. |
| `min_target_detection` | `0` | Do not remove condition-restricted targets in pooled Pando preprocessing. |
| `max_edges_per_target` | `Inf` | Preserve the complete structural candidate set. |

The Pando design is a structural hypothesis based on peak-to-gene domains and
motif evidence. Detection and regularized coefficient selection are handled by
RegCompass after this shared design is created.

### Multitask model

| Parameter | Default | Role |
|---|---:|---|
| `alpha` | `0.5` | Elastic-net balance between sparse selection and correlated-predictor stabilization. |
| `global_penalty_factor` | `1` | Explicit penalty factor for the shared backbone. |
| `deviation_penalty_factor` | `1` | Neutral explicit factor for zero-sum condition deviations. |
| `lambda_rule` | `"lambda.1se"` | Conservative lambda selected from leakage-resistant CV. |
| `nfolds` | `5` | Condition-stratified cell-level folds. |
| `n_bootstrap` | `100` | Full-size condition-stratified bootstrap replicates. |
| `min_selection_frequency` | `0.7` | Minimum bootstrap selection frequency. |
| `min_sign_stability` | `0.8` | Minimum conditional sign stability. |
| `min_bootstrap_success_fraction` | `0.8` | Minimum fraction of completed bootstrap fits. |
| `min_cv_rsq` | `0` | User floor; activation still requires strictly positive out-of-fold R-squared. |
| `candidate_screen_threshold` | `0` | Outcome-based full-data screening is disabled. |
| `max_edges_per_target` | `Inf` | No deterministic top-K truncation. |
| `min_detected_cells_per_condition` | `10` | Absolute observability floor. |
| `min_detection_fraction_per_condition` | `0.01` | Relative observability floor. |

## Shared structural and model edge universes

Pando returns a structural universe

\[
\mathcal U_m^{struct}
=\{e=(t,p,g):\text{motif and peak-to-gene rules support }e\}
\]

for cell type \(m\), using all conditions together.

RegCompass then defines, for condition \(c\),

\[
m_c=\min\left(
 n_c,
 \max\left(10,\left\lceil0.01n_c\right\rceil\right)
\right).
\]

For edge \(e=(t,p,g)\), let

\[
N^{TF\times peak}_{e,c}
=\sum_{u\in c}I(T_{t,u}>0\land A_{p,u}>0)
\]

and

\[
N^{target}_{e,c}
=\sum_{u\in c}I(Y_{g,u}>0).
\]

The model universe is

\[
\mathcal U_m^{model}
=\left\{e\in\mathcal U_m^{struct}:
\exists c,
N^{TF\times peak}_{e,c}\ge m_c,
N^{target}_{e,c}\ge m_c
\right\}.
\]

This rule is condition aware but outcome agnostic: it checks only whether the
actual interaction predictor and target are observable. It does not use their
correlation, coefficient, P value, or condition effect.

The structural and model universes receive separate MD5 fingerprints. Every
condition of the cell type uses exactly the same model-universe fingerprint.

## Multitask objective

For edge predictor

\[
x_{e,u}=T_{t,u}A_{p,u},
\]

condition-centred target and predictors are fitted with

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_c\delta_{e,c}=0.
\]

The weighted Gaussian elastic-net objective is

\[
\begin{aligned}
\min_{\beta,\gamma}\;&
\frac{1}{2}\sum_u w_u
\left(y_u-X_u\beta-X_uH_{c(u),\cdot}\gamma\right)^2\\
&+\lambda\alpha\sum_j v_j|b_j|
+\frac{\lambda(1-\alpha)}{2}\sum_jv_jb_j^2,
\end{aligned}
\]

where \(H=I-\mathbf 1\mathbf 1^T/C\), \(b=(\beta,\gamma)\), and

\[
w_u\propto\frac{1}{n_{c(u)}}.
\]

Each condition therefore contributes the same total regression loss.

### Why `alpha = 0.5`

TF–peak predictors are strongly correlated because nearby accessible peaks can
share motifs and TF expression is reused across several peaks. Pure lasso can
select one arbitrary member of a correlated set, whereas the ridge component of
an elastic net stabilizes correlated coefficients. `alpha = 0.5` is also the
current Pando regularized-model default and is retained as the canonical value.

### Why both penalty factors equal one

A fixed deviation factor of two was not a neutral default. In the symmetric
zero-sum parameterization, one biological deviation is represented across
multiple condition coordinates, and the deviation design columns have lower
norms than the global column. These properties already shrink deviations more
for a comparable fitted effect. Equal explicit factors avoid adding another
undocumented prior.

`deviation_penalty_factor > 1` remains valid as a prespecified sensitivity
analysis when the scientific prior explicitly favors a highly conserved shared
backbone. It must not be presented as universally optimal.

## Cross-validation

Five folds are assigned separately inside every condition. For each fold:

1. condition means for target and predictors are estimated on training cells;
2. training-derived means are applied to validation cells;
3. predictor scales are estimated on training cells;
4. the same scales are applied to validation cells;
5. lambda is evaluated on out-of-fold residuals.

The selected lambda uses the one-standard-error rule. A target model can define
active edges only when

\[
R^2_{CV}>0.
\]

Thus it must improve on the condition-centred zero predictor. A positive
`min_cv_rsq` adds a stronger user-specified floor.

With the canonical minimum of 100 cells per condition and five folds, each
condition contributes approximately 80 training and 20 validation cells per
fold. Lower overrides are allowed but generate a low-power warning when fewer
than 10 validation cells per condition are expected in each fold.

## Bootstrap reproducibility

Each bootstrap replicate samples \(n_c\) cells with replacement inside every
condition and refits at the full-data selected lambda.

For successful bootstrap set \(\mathcal B_s\),

\[
\widehat\Pi_{e,c}
=\frac{1}{B_s}\sum_{b\in\mathcal B_s}
I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon).
\]

The Monte Carlo standard error is

\[
SE(\widehat\Pi_{e,c})
=\sqrt{\frac{\widehat\Pi_{e,c}
(1-\widehat\Pi_{e,c})}{B_s}}.
\]

At \(\widehat\Pi=0.7\):

```text
B = 50  → SE ≈ 0.0648
B = 100 → SE ≈ 0.0458
```

Therefore the canonical default is 100. RegCompass records the Monte Carlo SE
and Wilson 95% interval for every condition edge.

Sign stability is

\[
\rho_{e,c}
=\left|
\frac{\sum_bI^{(b)}_{e,c}
\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_bI^{(b)}_{e,c}}
\right|.
\]

If \(q\) is the fraction of selected fits with the majority sign, then

\[
\rho=|2q-1|.
\]

Consequently, `min_sign_stability = 0.8` requires at least 90% agreement on one
sign among selected fits.

This procedure is described as full-size bootstrap reproducibility. It is not
the half-sample stability-selection algorithm and does not claim its formal
per-family error bound.

## Active-edge rule

An edge is active in condition \(c\) only when all conditions below hold:

\[
\begin{aligned}
& B_s/B\ge0.8,\\
& \widehat\Pi_{e,c}\ge0.7,\\
& \rho_{e,c}\ge0.8,\\
& |\widehat\theta_{e,c}|>\max(\tau_{effect},\varepsilon),\\
& R^2_{CV}>0,\\
& R^2_{CV}\ge\tau_{CV}\quad\text{when }\tau_{CV}>0.
\end{aligned}
\]

The Layer 1 projection coefficient is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}
\widehat\Pi_{e,c}\rho_{e,c}.
\]

## Supported sensitivity analyses

### Alternative peak-to-gene domain

Use GREAT as a separate sensitivity analysis:

```r
pando_args = list(
  pando_design_args = list(
    peak_to_gene_method = "GREAT",
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    min_tf_detection = 0,
    min_peak_detection = 0,
    min_target_detection = 0,
    max_edges_per_target = Inf
  )
)
```

Changing Signac to GREAT changes the structural edge universe, so Stage 1 and all
downstream stages must be rerun.

### Stronger shared-backbone prior

```r
multitask_args = list(
  deviation_penalty_factor = 2
)
```

This is a sensitivity prior, not the canonical default. Compare model-universe
fingerprints, CV R-squared, bootstrap completion, active target genes, and final
core reactions.

### Stricter predictive floor

```r
multitask_args = list(
  min_cv_rsq = 0.05
)
```

This changes active-edge membership but not the Pando structural candidate
universe.

## Unsupported shortcuts

The canonical model rejects:

- non-zero `candidate_screen_threshold`;
- finite `max_edges_per_target` in either Pando or multitask arguments;
- condition-specific Pando candidate dictionaries;
- using bootstrap edge frequencies as formal biological-replicate inference;
- interpreting deterministic candidate order as an evidence ranking.

## References used for the defaults

- Pando: TF-expression by peak-accessibility interaction models and current
  `infer_grn()` defaults (`alpha = 0.5`, Signac/GREAT domains, 100 kb upstream,
  0 downstream, and 1 Mb GREAT extension).
- Friedman, Hastie and Tibshirani: coordinate-descent elastic-net models.
- Meinshausen and Bühlmann: stability-selection frequency thresholds and the
  practical use of 100 repeated selections. RegCompass uses a different
  full-size bootstrap and therefore reports reproducibility without claiming
  the half-sample error-control theorem.
- SCENIC gene filtering: detection in at least 1% of cells as a soft noise
  filter. RegCompass adds a 10-cell floor and applies it condition-wise to the
  actual TF×peak predictor and target.
