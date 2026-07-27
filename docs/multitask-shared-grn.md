# Shared-background multitask GRN mathematics and object contracts

RegCompassR 1.8.9 separates the structural regulatory hypothesis, the observable shared model dictionary, and the condition-specific fitted network.

## 1. Shared structural background

For each cell type `m`, all conditions are combined only to construct one Pando structural dictionary:

\[
\mathcal U_m^{struct}=\{e=(t,p,g)\},
\]

where `t` is a TF, `p` is an exact measured ATAC feature overlapping a Pando regulatory region containing a motif for `t`, and `g` is a GEM GPR target gene linked to that region.

`Pando::prepare_grn_design()` returns one validated version-2 `PandoGRNDesign` per cell type. Candidate edges are deduplicated by `(TF, ATAC feature, target)`. Multiple regulatory regions mapping to one measured peak remain in `supporting_regions` but do not duplicate a design-matrix column. The structural design MD5 fingerprint is stored as `edge_universe_id`.

Canonical Pando pooled detection thresholds are zero. No pooled fitted coefficient, pooled p-value, or pooled target-correlation screen is required for candidate entry, because true opposite condition effects can cancel and condition-restricted regulators can be diluted in pooled data.

## 2. Condition-aware observable model universe

For condition `c` with `n_c` cells, define

\[
m_c=\min\left(
 n_c,
 \max\left(10,\left\lceil0.01n_c\right\rceil\right)
\right).
\]

For edge \(e=(t,p,g)\),

\[
N^{TF\times peak}_{e,c}
=\sum_{u\in c}I(T_{t,u}>0\land A_{p,u}>0)
\]

and

\[
N^{target}_{e,c}
=\sum_{u\in c}I(Y_{g,u}>0).
\]

The shared model universe is

\[
\mathcal U_m^{model}
=\left\{e\in\mathcal U_m^{struct}:
\exists c,
N^{TF\times peak}_{e,c}\ge m_c,
N^{target}_{e,c}\ge m_c
\right\}.
\]

TF and peak must be non-zero in the same cells because the fitted predictor is their product. The filter uses detection only and therefore does not leak target correlation, coefficient magnitude, or condition effect into candidate selection.

The model universe receives a second MD5 fingerprint, `model_edge_universe_id`. Every condition of the cell type uses the same model fingerprint. The full structural candidate table retains `model_observable`, `n_observable_conditions`, and `observable_conditions` so filtering remains auditable.

Finite `max_edges_per_target` values are rejected in canonical mode. Pando candidate order is deterministic but is not a statistical evidence ranking.

## 3. Predictor, condition centring, and shared scale

For cell `u` and edge `e=(t,p,g)`:

\[
x_{e,u}=T_{t,u}A_{p,u},
\]

where `T` is normalized TF RNA and `A` is cell-type-shared TF-IDF peak accessibility.

Target expression and predictors are centred only within condition:

\[
y^\circ_{g,u}=y_{g,u}-\bar y_{g,c(u)},
\qquad
x^\circ_{e,u}=x_{e,u}-\bar x_{e,c(u)}.
\]

The shared edge scale is

\[
s_{e,m}=
\sqrt{
\frac1C\sum_{c=1}^{C}
\frac1{n_c}\sum_{u\in c}(x^\circ_{e,u})^2
},
\qquad
\widetilde x_{e,u}=x^\circ_{e,u}/s_{e,m}.
\]

The same scale is used for every condition. The canonical public workflow does not accept or interpret a biological-sample metadata column.

## 4. Global backbone and symmetric condition deviations

RegCompass estimates

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_{c=1}^{C}\delta_{e,c}=0.
\]

The symmetric centring matrix is

\[
H=I_C-\frac1C\mathbf 1\mathbf 1^T,
\]

and the design contains

\[
[\widetilde X,\;\widetilde XH_1,\ldots,\widetilde XH_C].
\]

The elastic-net objective is

\[
\begin{aligned}
\min_{\beta,\gamma}\quad
&\frac12\sum_u w_u
\left[
 y^\circ_u-
 \widetilde X_u\beta-
 \sum_{k=1}^{C}H_{c(u),k}\widetilde X_u\gamma_k
\right]^2\\
&+\lambda\alpha
\left[p_G\|\beta\|_1+p_D\sum_k\|\gamma_k\|_1\right]\\
&+\frac{\lambda(1-\alpha)}2
\left[p_G\|\beta\|_2^2+p_D\sum_k\|\gamma_k\|_2^2\right].
\end{aligned}
\]

The deviation is

\[
\delta=H\gamma,
\]

so

\[
\frac1C\sum_c\theta_{e,c}=\beta_e.
\]

Canonical parameters are

\[
\alpha=0.5,\qquad p_G=1,\qquad p_D=1.
\]

`alpha = 0.5` retains a lasso component for sparse edge selection and a ridge component that stabilizes correlated TF–peak predictors and the rank-deficient symmetric deviation representation. Global and deviation coordinates receive equal explicit penalty factors by default. The zero-sum deviation block already has multiplicity and lower design-column norms; an additional fixed factor of two is therefore a stronger prior rather than a neutral default. `p_D > 1` remains an explicit sensitivity analysis favoring a conserved shared backbone.

Observation weights are

\[
w_u\propto\frac1{n_{c(u)}},
\]

so every condition contributes equal total loss regardless of cell count.

## 5. Leakage-resistant cross-validation

The canonical model uses five folds stratified within every condition. For each fold:

1. condition means for the target and predictors are estimated from training cells;
2. training means are applied to validation cells;
3. predictor scales are estimated from training cells;
4. training scales are applied to validation cells;
5. every lambda is evaluated on out-of-fold residuals.

The selected lambda uses `lambda.1se`. The reported

\[
R^2_{CV}=1-
\frac{\sum_uw_u(y_u^\circ-\widehat y_u^\circ)^2}
     {\sum_uw_u(y_u^\circ-\bar y_w^\circ)^2}
\]

is an internal predictive reliability statistic for the joint target model. An edge can be active only when

\[
R^2_{CV}>0.
\]

Thus the model must improve on the condition-centred null predictor. A positive `min_cv_rsq` adds a stronger user floor. `cv_rsq` is not a biological-replicate significance test.

The canonical minimum is 100 cells per condition and cell type. With five folds this yields about 80 training and 20 validation cells per condition in each fold. A lower override is allowed but warns when fewer than 10 validation cells per condition are expected.

## 6. Condition-stratified bootstrap reproducibility

For every bootstrap replicate, condition `c` is sampled with replacement using exactly its original cell count `n_c`. The bootstrap target and predictor matrices are re-centred within the resampled condition:

\[
y^{\circ(b)}_{g,u}
=y^{(b)}_{g,u}-\bar y^{(b)}_{g,c(u)},
\]

\[
x^{\circ(b)}_{e,u}
=x^{(b)}_{e,u}-\bar x^{(b)}_{e,c(u)}.
\]

Re-centring is required because duplicated bootstrap rows do not preserve the zero means of the full-data residual matrix. Each bootstrap uses the full-data shared scale `s_e` and selected lambda, so coefficients remain comparable across replicates.

For successful bootstrap set \(\mathcal B_s\), selection frequency is

\[
\widehat\Pi_{e,c}
=\frac1{|\mathcal B_s|}
\sum_{b\in\mathcal B_s}
I(|\theta^{(b)}_{e,c}|>\varepsilon).
\]

Its Monte Carlo standard error is

\[
SE(\widehat\Pi_{e,c})
=\sqrt{\frac{\widehat\Pi_{e,c}(1-\widehat\Pi_{e,c})}
{|\mathcal B_s|}}.
\]

RegCompass also reports a Wilson 95% interval. At \(\widehat\Pi=0.7\), the Monte Carlo standard error is approximately 0.0648 for 50 replicates and 0.0458 for 100 replicates. The canonical default is therefore `n_bootstrap = 100L`.

Conditional sign stability is

\[
\rho_{e,c}
=\left|
\frac{\sum_b I^{(b)}_{e,c}
\operatorname{sign}(\theta^{(b)}_{e,c})}
{\sum_b I^{(b)}_{e,c}}
\right|.
\]

When `q` is the fraction of selected fits with the majority sign,

\[
\rho=|2q-1|.
\]

Therefore `min_sign_stability = 0.8` requires at least 90% agreement on one sign among selected fits.

The reliability-weighted coefficient used for Layer 1 is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}
\widehat\Pi_{e,c}\rho_{e,c}.
\]

The full-size bootstrap estimates cell-resampling reproducibility. It is not the half-sample stability-selection procedure and does not claim its formal per-family error bound.

## 7. Active-edge rule

An edge is active only when all conditions hold:

\[
\begin{aligned}
&|\mathcal B_s|/B\ge0.8,\\
&\widehat\Pi_{e,c}\ge0.7,\\
&\rho_{e,c}\ge0.8,\\
&|\widehat\theta_{e,c}|>\max(\tau_{effect},\varepsilon),\\
&R^2_{CV}>0,\\
&R^2_{CV}\ge\tau_{CV}\quad\text{when }\tau_{CV}>0.
\end{aligned}
\]

Multitask mode does not fabricate a classical adjusted p-value: `padj` is `NA`, and `evidence_type` is `multitask_bootstrap_stability_selected`.

## 8. Condition genes and complete-GPR cores

The condition sub-GRN target set is

\[
G_{m,c}=\{g:\exists e=(t,p,g)\text{ active in }c\}.
\]

Both positive and negative active edges establish membership because both indicate condition-specific regulation.

For reaction `r` with alternative GPR branches `B_{r,k}`:

\[
Core_{r,m,c}=1
\iff
\exists k:B_{r,k}\subseteq G_{m,c}.
\]

This preserves the distinction between required complex subunits and alternative isoenzymes.

## 9. ATAC-only regulatory projection

The GRN is fitted with RNA and ATAC, but the downstream modifier reprojects only metacell ATAC. For edge `(t,p,g)`:

\[
\omega_{tpg,c}
=\widetilde\theta_{tpg,c}
\frac{\bar T_{t,m}}{s_{tpg,m}},
\]

where `bar T` is the equal-condition TF reference from Stage 1. This is the local ATAC derivative of the standardized TF-by-ATAC predictor evaluated at one shared TF reference.

TFs sharing the same measured peak and target are signed-summed:

\[
\psi_{p,g,c}=\sum_t\omega_{tpg,c}.
\]

Let `D_{p,u}` be the bounded robust ATAC deviation relative to all metacells of the same cell type. One denominator is shared across conditions:

\[
L_{g,m}=\max_c\sum_p|\psi_{p,g,c}|.
\]

The modifier is

\[
R_{g,u,c}
=q_{g,m}
\operatorname{clip}
\left(
\frac{\sum_p\psi_{p,g,c}D_{p,u}}{L_{g,m}+\varepsilon},
-1,1
\right),
\]

with

\[
q_{g,m}=\sqrt{\operatorname{clip}(R^2_{g,m,CV},0,1)}.
\]

No active edge in one condition gives `R = 0`, producing exact RNA-only fallback:

\[
C^{MO}_{g,u,c}
=\frac{C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}
{1-C^{RNA}_{g,u,c}+C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}.
\]

## 10. Shared final metabolic structure

Condition core sets may differ:

\[
Core_{m,A}\neq Core_{m,B}.
\]

Stage 3 merges all biological reaction catalogues. For each medium, Stage 5 constructs one union GEM and performs one global FASTCORE completion:

\[
S_A=S_B=S^{union},
\qquad
lb_A=lb_B,
\qquad
ub_A=ub_B.
\]

Only evidence-derived reaction penalties differ by metacell. This separates regulatory/metabolic evidence differences from structural model differences.

## 11. Stage object contract

Stage 1 exposes:

```text
tf_peak_gene_candidates
tf_peak_gene_global
tf_peak_gene_condition_all
tf_peak_gene_significant
condition_target_genes
target_model_diagnostics
stability_diagnostics
group_status
celltype_fit_status
```

`tf_peak_gene_candidates` retains the complete Pando structural universe and identifies edges entering the model with `model_observable` and `model_edge_universe_id`.

`tf_peak_gene_significant` retains the columns needed by Stage 3 (`group_id`, condition, cell type, `tf`, `region`, `target`, `estimate`, `rsq`) and adds global, deviation, effective, CV, bootstrap, Monte Carlo uncertainty, and stability fields.

Stage 3 exposes condition supported genes, complete-GPR cores, and biological reaction membership using `group_id` as the actual `condition × cell type` analysis unit. The merged catalogue is not itself a GEM. Stage 5 alone constructs the shared medium-specific union GEM.

See [Pando and multitask GRN parameter policy](grn-parameter-policy.md) for the canonical defaults and supported sensitivity analyses.
