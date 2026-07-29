# Pando shared-design, condition-comparable GRNs

> Each broad cell type is trained, validated, and refitted independently. If `cell_type` is supplied, only that label or those labels are processed. OOF folds are condition-stratified cells from the same fitted type; no cells from another type enter training or validation. Biological sample metadata and sample count are not inputs or gates.

RegCompass Stage 1 calls `Pando::initiate_grn()`, `Pando::find_motifs()`, and
`Pando::infer_condition_grn()` once on the normalized paired-cell multiome
object. Pando returns one versioned `ConditionGRNFit` per broad cell type.

## Condition-sparse multitask common-metric model

The canonical configuration is:

```r
pando_infer_args = list(
  candidate_screen = "motif_domain",
  condition_mix = 0.5,
  condition_weight = "equal",
  reference_condition = "Control",
  scale = TRUE
)
```

For a target gene, let `X_c` and `y_c` denote the predictor matrix and response
for condition `c`, and let `B` contain one coefficient column per condition.
For fixed `lambda`, `alpha`, and `rho = condition_mix`, Pando minimizes the
coupled sparse-group multitask objective

\[
\sum_c \frac{1}{2n_c}\lVert y_c-a_c-X_c\beta_c\rVert_2^2
+\frac{\lambda(1-\alpha)}{2}\lVert B\rVert_F^2
+\lambda\alpha\left[
(1-\rho)\sum_e\lVert B_{e,\cdot}\rVert_2
+\rho\sum_{e,c}|\beta_{e,c}|
\right].
\]

The condition columns are therefore **not independent elastic-net fits**. The
row-wise group penalty couples selection of the same TF–peak–target edge across
conditions, whereas the element-wise L1 term permits condition-specific zeros.
After sparse support selection, Pando performs a support-constrained ridge refit
using a pooled common predictor metric and an estimated shared baseline.
Condition coefficients may still differ in magnitude or sign, but they are
estimated under this explicit cross-condition coupling.

The following quantities are shared for direct comparison:

- one TF–peak–target coordinate system;
- one equal-condition, within-condition-variance transform for each interaction predictor, learned inside each outer training fold;
- one equal-condition target transform;
- one lambda path and one selected lambda per target, selected from mean condition validation loss;
- one pooled common predictor metric for the support-constrained refit;
- explicit edge-by-condition estimability, selected-support, and active masks.

## Candidate-edge policy

The fitted predictor is

\[
x_{e,i}=RNA_{TF(e),i}\times ATAC_{peak(e),i}.
\]

Marginal TF–target and peak–target correlations can be small even when this
interaction predicts the target. RegCompass therefore defaults to
`candidate_screen = "motif_domain"`: motif/domain candidates are retained
without response-dependent screening, and sparse-group regularization performs
coefficient selection.

One optional Pando mode remains available:

- `pooled_within_condition`: apply pooled within-condition marginal association screening before model fitting.

This is a sensitivity/performance mode, not the canonical interaction-safe
default. Because the response is used in candidate screening outside the nested
OOF loop, Pando marks its projection as ineligible for primary penalty
construction.

## Equal-condition coefficient units

For `K` fitted conditions, Pando assigns each condition weight `1 / K`,
independent of its cell count. For each edge predictor it computes

\[
\mu_e=\frac{1}{K}\sum_c \bar{x}_{e,c},
\]

\[
s_e^2=\frac{1}{K}\sum_c
\left[\frac{1}{n_c}\sum_{i\in c}(x_{e,i}-\bar{x}_{e,c})^2\right],
\]

and applies

\[
z_{e,i}=\frac{x_{e,i}-\mu_e}{s_e}.
\]

Thus the scale is the equal-condition average of within-condition population
variances; it does not include between-condition mean separation. The target
response uses the analogous equal-condition center and within-condition
variance scale. Each outer OOF fold learns these quantities from its training
cells only and applies the stored transform to held-out cells.

Scaling the final interaction is not equivalent to separately scaling RNA and
ATAC before multiplication.

The `ConditionGRNFit` v5 contract retains, among other fields:

```text
edge_table
beta_condition
beta_shared
delta_condition
contrast
estimability_mask
support_mask
active_mask
comparison_mask
predictor_transform
response_transform
intercept
target_fit
condition_rsq_oof
target_rsq_oof_pooled
projection_common_oof
projection_condition_full_oof
projection_global_common_oof
```

## Reference contrasts and comparison support

For reference condition `r`, Pando retains the interpretation-layer contrast

\[
\Delta\beta_{e,c}=\beta_{e,c}-\beta_{e,r}.
\]

A zero or unavailable coefficient has distinct semantics:

1. an estimable coefficient selected or refitted as zero remains numeric zero;
2. an edge not estimable in a condition is represented as `NA` in `beta_condition`.

Pando retains the reference-comparison relationship with

\[
comparison\_mask_{e,c}=
 estimability_{e,c}\land estimability_{e,r}.
\]

RegCompass validates the dimensions, dimnames, logical type, reference
condition, and exact relationship to `estimability_mask`. It does not silently
reconstruct missing comparison support from an older Pando object.

The complete condition-effect table retains all rows and adds
`comparable_to_reference`. Rows with `TRUE` may enter the active
**reference-effect table**. This table is an interpretation artifact: the
canonical Stage 4 primary penalty does not project `Delta beta` and does not
consume the active reference-effect table.

## Stage 3 use

Stage 3 uses active absolute condition coefficients `beta[e, c]` to identify
supported metabolic target genes and complete-GPR core reactions. It does not
require a reference contrast.

## Stage 4 primary penalty projection

The canonical Stage 4 path requires:

```text
projection_component = "condition"
origin = "oof"
support_policy = "pairwise_common"   # two-condition comparison
                 or "global_common"  # multi-condition comparison
```

For held-out cell `i` in condition `c`, outer fold `k`, and target gene `g`, the
stored common-support Pando projection is

\[
G^{OOF}_{i,g,c}=\sum_{e:\,target(e)=g,\,m_e=1}
\beta^{(-k)}_{e,c}z^{(-k)}_{i,e},
\]

where `m_e` is the pairwise-common or global-common estimability mask. These are
absolute condition coefficients, not reference contrasts. Pando computes the
TF-RNA × peak-ATAC product, applies the training-fold transform, and projects
the held-out cell before RegCompass receives the result.

RegCompass aggregates the already-computed single-cell scores by arithmetic
mean over the SuperCell membership:

\[
G^{OOF}_{u,g,c}=\frac{1}{|\mathcal{M}_u|}
\sum_{i\in\mathcal{M}_u}G^{OOF}_{i,g,c}.
\]

It does not recompute TF × ATAC from metacell means, refit coefficients,
recenter or rescale by condition, or convert unavailable projections to zero.

The `condition_estimable` projection, which can retain condition-specific
estimable edges, is produced only as a diagnostic comparator and is explicitly
ineligible for the primary cross-condition penalty.

### Pooled OOF reliability

For each target gene, Pando stores

\[
R^2_{OOF,pooled,g}=1-
\frac{\sum_c\sum_{i\in c}(y_{i,g}-\hat y^{OOF}_{i,g})^2}
{\sum_c\sum_{i\in c}(y_{i,g}-\bar y_{c,g})^2}.
\]

RegCompass broadcasts the condition-pooled reliability weight

\[
q_g=\sqrt{clamp(R^2_{OOF,pooled,g},0,1)}
\]

to every metacell covered by that broad-cell-type fit. A non-finite pooled OOF
`R^2`, unavailable predictive OOF status, or incomplete OOF cell coverage makes
the reliability unavailable.

### Broad-cell-type pooled calibration scale

For each target gene and broad cell type, RegCompass computes one calibration
scale from the finite common-support projections of **all conditions** in that
cell type:

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},
10^{-6}
\right).
\]

If fewer than two finite values are available or the projection has no range,
the fallback scale is `1`.

The signed bounded modifier is

\[
R_{g,c,u}=q_g\tanh\left(\frac{G^{OOF}_{u,g,c}}{\sigma_{g,t}}\right).
\]

Therefore the current Stage 4 implementation uses:

- absolute condition coefficients, not `Delta beta`;
- pairwise-common or global-common estimable edges for the primary penalty;
- condition-pooled OOF reliability;
- one broad-cell-type pooled robust projection scale;
- no absolute-sum coefficient normalization.

Because `sigma[g, t]` is estimated from the projection distribution, this
modifier is **not 1-homogeneous** in the coefficients. Multiplying all
projections by the same positive constant also multiplies the pooled scale, so
`G / sigma` generally remains unchanged rather than scaling with the
coefficients.

The modifier is integrated with bounded RNA support on the log-odds scale:

\[
C_{multiome}=\frac{C_{RNA}\,2^{\alpha R}}
{1-C_{RNA}+C_{RNA}\,2^{\alpha R}}.
\]

RNA support values exactly equal to `0` or `1` remain absorbing values.
RegCompass then applies the shared GPR rules, expression-linked reaction
penalty, shared GEM, medium, reaction bounds, and directional two-step LP.

## Pando bridge ownership

RegCompass forwards three argument bundles:

```text
pando_initiate_args → initiate_grn
pando_motif_args    → find_motifs
pando_infer_args    → infer_condition_grn
```

It manages the object, assays, motif object, genome, metadata columns, complete
GEM target-gene set, network name, minimum condition size, error policy, and
`BPPARAM`. Nested attempts to override these fields stop before Pando is called.
Canonical Stage 1 also rejects Pando aggregation columns because Stage 2 owns
metacell construction.

## Persisted artifacts

Stage 1 writes:

```text
pando_group_status.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_condition_active.tsv.gz
pando_tf_peak_gene_condition_effect_all.tsv.gz
pando_tf_peak_gene_condition_effect_active.tsv.gz
pando_tf_peak_gene_universal.tsv.gz
pando_condition_network_index.tsv.gz
pando_condition_fit_diagnostics.tsv.gz
pando_edge_predictor_transforms.tsv.gz
pando_condition_grn_fits.rds
pando_objects/condition_grn_fit_v5.rds
```

The in-memory audit surface is:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_condition_effect_all
step1$grn_result$normalization_policy
step1$params$pando_parallel
```

## Genome-build safety

Human analyses may use the bundled hg38 phastCons plus SCREEN ccRE union. Mouse
analyses must supply a species- and build-matched `GRanges` through
`pando_initiate_args$regions`. RegCompass stops rather than applying hg38
regulatory regions to mouse ATAC coordinates.
