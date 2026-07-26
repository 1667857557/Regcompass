# Shared-background multitask GRN mathematics and object contracts

RegCompassR 1.8.8 separates the structural regulatory hypothesis from the condition-specific fitted network.

## 1. Shared structural background

For each cell type `m`, all conditions are combined only to construct one candidate dictionary:

\[
\mathcal U_m=\{e=(t,p,g)\},
\]

where `t` is a TF, `p` is an exact measured ATAC feature overlapping a Pando regulatory region containing a motif for `t`, and `g` is a GEM GPR target gene linked to that region.

`Pando::prepare_grn_design()` returns one validated version-2 `PandoGRNDesign` per cell type. Candidate edges are deduplicated by `(TF, ATAC feature, target)`. Multiple regulatory regions mapping to one measured peak remain in `supporting_regions` but do not duplicate a design-matrix column. The design MD5 fingerprint is validated before fitting.

No pooled fitted coefficient or pooled p-value is required for candidate entry, because true opposite condition effects can cancel under pooled screening.

## 2. Predictor, condition centring and shared scale

For cell `u` and edge `e=(t,p,g)`:

\[
x_{e,u}=T_{t,u}A_{p,u},
\]

where `T` is normalised TF RNA and `A` is cell-type-shared TF-IDF peak accessibility.

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

## 3. Global backbone and symmetric condition deviations

RegCompass estimates

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_{c=1}^{C}\delta_{e,c}=0.
\]

The symmetric centring matrix is

\[
H=I_C-\frac1C\mathbf 1\mathbf 1^T.
\]

The design contains

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

The condition deviation is

\[
\delta=H\gamma.
\]

Therefore

\[
\frac1C\sum_c\theta_{e,c}=\beta_e.
\]

Because the symmetric `C`-column deviation representation is rank deficient, RegCompass requires `alpha < 1`; the non-zero ridge component makes the regularised solution unique.

Observation weights are

\[
w_u\propto\frac1{n_{c(u)}},
\]

so every condition contributes equal total loss regardless of cell count.

## 4. Cross-validation

`lambda.min` or `lambda.1se` is selected by folds stratified within condition at cell level. The reported `cv_rsq` is an internal out-of-fold reliability statistic for the joint target model. It is not a biological-replicate significance test.

## 5. Condition-stratified bootstrap stability

For every bootstrap replicate, condition `c` is sampled with replacement using exactly its original cell count `n_c`. The bootstrap target and predictor matrices are then re-centred within the resampled condition:

\[
y^{\circ(b)}_{g,u}
=y^{(b)}_{g,u}-\bar y^{(b)}_{g,c(u)},
\]

\[
x^{\circ(b)}_{e,u}
=x^{(b)}_{e,u}-\bar x^{(b)}_{e,c(u)}.
\]

Re-centring is required because duplicated bootstrap rows do not preserve the zero means of the full-data residual matrix. Each bootstrap uses the full-data shared scale `s_e` and the full-data selected `lambda`, so coefficients remain comparable across replicates.

Selection frequency is

\[
\Pi_{e,c}
=\frac1B\sum_b I(\theta^{(b)}_{e,c}\neq0).
\]

Conditional sign stability is

\[
\rho_{e,c}
=
\left|
\frac{\sum_b I(\theta^{(b)}_{e,c}\neq0)
\operatorname{sign}(\theta^{(b)}_{e,c})}
{\sum_b I(\theta^{(b)}_{e,c}\neq0)}
\right|.
\]

The reliability-weighted coefficient used for Layer 1 is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

An active edge must separately satisfy the configured selection-frequency, sign-stability, absolute-effect, and `cv_rsq` thresholds. Multitask mode does not fabricate a classical adjusted p-value: `padj` is `NA`, and `evidence_type` is `multitask_bootstrap_stability_selected`.

The default is `n_bootstrap = 50L`; `100` or more is recommended for a final analysis when computation permits.

## 6. Condition genes and complete-GPR cores

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

## 7. ATAC-only regulatory projection

The GRN is fitted with RNA and ATAC, but the downstream modifier reprojects only metacell ATAC. For edge `(t,p,g)`:

\[
\omega_{tpg,c}
=
\widetilde\theta_{tpg,c}
\frac{\bar T_{t,m}}{s_{tpg,m}},
\]

where `bar T` is the equal-condition TF reference from Stage 1. This is the local ATAC derivative of the standardised TF-by-ATAC predictor evaluated at one shared TF reference.

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
=
q_{g,m}
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
=
\frac{C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}
{1-C^{RNA}_{g,u,c}+C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}.
\]

## 8. Shared final metabolic structure

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

## 9. Stage object contract

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

`tf_peak_gene_significant` retains the columns needed by Stage 3 (`group_id`, condition, cell type, `tf`, `region`, `target`, `estimate`, `rsq`) and adds global, deviation, effective, bootstrap, and stability fields.

Stage 3 exposes condition supported genes, complete-GPR cores, and biological reaction membership using `group_id` as the actual `condition × cell type` analysis unit. The merged catalogue retains candidate-universe IDs and is not itself a GEM. Stage 5 alone constructs the shared medium-specific union GEM.
