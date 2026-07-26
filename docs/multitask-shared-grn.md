# Shared-background multitask GRN mathematics and object contracts

RegCompassR 1.8.8 separates the structural regulatory hypothesis from the condition-specific fitted network.

## 1. Shared structural background

For each cell type `m`, all conditions are pooled only for construction of the candidate dictionary:

\[
\mathcal U_m=\{e=(t,p,g)\},
\]

where `t` is a TF, `p` is a measured ATAC peak overlapping a Pando regulatory region containing a motif for `t`, and `g` is a GEM GPR target gene linked to that region.

`Pando::prepare_grn_design()` returns one `PandoGRNDesign` object per cell type. Candidate edges are deduplicated by the exact predictor identity `(TF, ATAC feature, target)`. Multiple Pando regulatory regions mapping to the same measured peak are retained as provenance in `supporting_regions`, but they do not duplicate a design-matrix column.

The candidate universe is condition agnostic. No pooled coefficient or pooled p-value is required for entry, because opposite condition effects could cancel under pooled screening.

## 2. Predictor and residualisation

For cell `u` and edge `e=(t,p,g)`:

\[
x_{e,u}=T_{t,u}A_{p,u},
\]

where `T` is normalised TF RNA and `A` is cell-type-shared TF-IDF peak accessibility.

When no biological sample column is supplied, target expression and predictors are centred within condition:

\[
y^\circ_{g,u}=y_{g,u}-\bar y_{g,c(u)},
\qquad
x^\circ_{e,u}=x_{e,u}-\bar x_{e,c(u)}.
\]

When `sample_col` is supplied, centring uses `condition × sample`, preventing sample-level baseline shifts from being interpreted as within-sample regulatory coupling.

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

The same scale is used for every condition.

## 3. Global backbone and symmetric condition deviations

RegCompass estimates

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_{c=1}^{C}\delta_{e,c}=0.
\]

The implementation uses the symmetric centring matrix

\[
H=I_C-\frac1C\mathbf 1\mathbf 1^T.
\]

For condition indicator row `z_u`, the design contains:

\[
[\widetilde X,\; \widetilde X H_1,\ldots,\widetilde X H_C].
\]

All condition deviation blocks receive the same penalty. The elastic-net objective is

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

The condition deviation is decoded as

\[
\delta=H\gamma.
\]

Therefore the deviations sum to zero and the arithmetic mean of the condition coefficients equals the global coefficient:

\[
\frac1C\sum_c\theta_{e,c}=\beta_e.
\]

Because the centred `C`-column representation is rank deficient, RegCompass requires `alpha < 1`. The non-zero ridge component makes the symmetric coefficient solution unique.

Observation weights are

\[
w_u\propto\frac1{n_{c(u)}},
\]

so every condition contributes the same total loss weight regardless of cell count.

## 4. Cross-validation and stability

`lambda.min` or `lambda.1se` is selected by condition-stratified folds. When each condition has at least three biological samples, complete samples are assigned to folds; otherwise folds are stratified at cell level.

The reported `cv_rsq` is an internal out-of-fold reliability statistic for the joint target model. It is not an independent biological-replicate significance test.

For repeated stratified subsampling:

\[
\Pi_{e,c}
=\frac1B\sum_b I(\theta^{(b)}_{e,c}\neq0),
\]

\[
\rho_{e,c}
=
\left|
\frac{\sum_b I(\theta^{(b)}_{e,c}\neq0)
\operatorname{sign}(\theta^{(b)}_{e,c})}
{\sum_b I(\theta^{(b)}_{e,c}\neq0)}
\right|.
\]

The stable coefficient is

\[
\widetilde\theta_{e,c}
=\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

An active edge must satisfy the configured selection-frequency, sign-stability, absolute-effect and model-reliability thresholds. Multitask mode does not create a classical adjusted p-value; `padj` is intentionally `NA` and `evidence_type` records stability selection.

## 5. Condition genes and reactions

The condition sub-GRN target set is

\[
G_{m,c}=\{g:\exists e=(t,p,g)\text{ active in }c\}.
\]

Both positive and negative active edges establish membership because both indicate condition-specific regulation.

For reaction `r` with GPR branches `B_{r,k}`:

\[
Core_{r,m,c}=1
\iff
\exists k:B_{r,k}\subseteq G_{m,c}.
\]

This retains the distinction between required complex subunits and alternative isoenzymes.

## 6. ATAC-only regulatory projection

The GRN is fitted with RNA and ATAC, but the downstream modifier reprojects only metacell ATAC. For each fitted edge:

\[
\omega_{tpg,c}
=
\widetilde\theta_{tpg,c}
\frac{\bar T_{t,m}}{s_{tpg,m}},
\]

where `bar T` is the equal-condition TF reference computed during Stage 1.

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

A condition with no active edge has `R = 0` and therefore exact RNA-only fallback under the existing log-odds integration:

\[
C^{MO}_{g,u,c}
=
\frac{C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}
{1-C^{RNA}_{g,u,c}+C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}.
\]

## 7. Shared final metabolic structure

Condition core sets may differ:

\[
Core_{m,A}\neq Core_{m,B}.
\]

Stage 3 merges their biological reaction catalogues. For each medium, Stage 5 constructs one union GEM and performs one global FASTCORE completion:

\[
S_A=S_B=S^{union},
\qquad
lb_A=lb_B,
\qquad
ub_A=ub_B.
\]

Only reaction penalties differ by metacell. This separates biological evidence differences from structural model differences.

## 8. Stage object contract

Stage 1 exposes:

```text
tf_peak_gene_candidates
tf_peak_gene_global
tf_peak_gene_condition_all
tf_peak_gene_significant
condition_target_genes
target_model_diagnostics
stability_diagnostics
sample_status
celltype_fit_status
```

`tf_peak_gene_significant` retains the legacy columns needed by Stage 3 (`group_id`, condition, cell type, `tf`, `region`, `target`, `estimate`, `rsq`) while adding global, deviation, effective and stability fields.

Stage 3 exposes condition supported genes, complete-GPR cores and biological reaction membership. The merged catalogue retains candidate-universe IDs and is not itself a GEM. Stage 5 is the only stage that constructs the shared medium-specific union GEM.
