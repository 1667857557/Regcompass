# Pando shared-design, condition-comparable GRNs

> Each broad cell type is trained, validated, and refitted independently. If `cell_type` is supplied, only that label or those labels are processed. OOF folds are condition-stratified cells from the same fitted type; no cells from another type enter training or validation. Biological sample metadata and sample count are not inputs or gates.

RegCompass Stage 1 calls `Pando::initiate_grn()`, `Pando::find_motifs()`, and
`Pando::infer_condition_grn()` once on the normalized paired-cell multiome
object. Pando returns one versioned `ConditionGRNFit` per cell type.

## Condition-sparse common-metric model

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

For target gene `g`, condition `c`, edge matrix `X_c`, and target-specific shared
lambda, Pando solves an ordinary elastic-net problem for each condition:

\[
\frac{1}{2n_c}\lVert y_c-a_c-X_c\beta_c\rVert_2^2+
\lambda\left[\eta\lVert\beta_c\rVert_1+
\frac{1-\eta}{2}\lVert\beta_c\rVert_2^2\right].
\]

Condition coefficient columns are separable at fixed lambda. Cross-condition
sharing is restricted to quantities required for direct comparison:

- one TF–peak–target coordinate system;
- one equal-condition, within-condition-variance transform of each interaction predictor, learned inside each outer fold;
- one pooled target transform;
- one lambda path and selected lambda per target;
- explicit edge-by-condition estimability metadata.

## Candidate-edge policy

The fitted predictor is

\[
x_{e,u}=RNA_{TF(e),u}\times ATAC_{peak(e),u}.
\]

Marginal TF-target and peak-target correlations can be small even when this
interaction predicts the target. RegCompass therefore defaults to
`candidate_screen = "motif_domain"`: structural motif/domain candidates are
retained, and elastic-net regularization performs coefficient selection.

One optional Pando mode remains available:

- `pooled_within_condition`: remove condition means, then apply marginal TF-target and peak-target screening before model fitting.

This is a sensitivity/performance mode, not the canonical interaction-safe default. Because the response is used before cross-validation, its projection is ineligible for penalty construction.

## Common coefficient units

Pando computes one center `mu_e` and scale `s_e` over all cells of the cell type:

\[
z_{e,u}=\frac{x_{e,u}-\mu_e}{s_e}.
\]

The target is also centered and scaled once over the pooled cell type. Scaling
the final interaction is not equivalent to separately scaling RNA and ATAC
before multiplication.

The `ConditionGRNFit` retains:

```text
edge_table
beta
contrast
eligibility_mask
comparison_mask
predictor_transform
response_transform
intercept
target_fit
target_rsq
lambda_path
cv
```

## Reference contrasts and comparison support

For reference condition `r`:

\[
\Delta\beta_{e,c}=\beta_{e,c}-\beta_{e,r}.
\]

A zero coefficient can have two distinct meanings:

1. an eligible coefficient estimated as zero by elastic net;
2. a coefficient fixed to zero because the edge is not estimable in that
   condition.

Pando 1.5.0 retains this interpretation-layer distinction with:

\[
comparison\_mask_{e,c}=
 eligibility_{e,c}\land eligibility_{e,r}.
\]

RegCompass requires this explicit field. It validates the matrix dimensions,
dimnames, logical type, reference condition, and exact relationship to
`eligibility_mask`. It does not silently reconstruct missing comparison support
from an older Pando object.

The complete condition-effect table retains all rows and adds
`comparable_to_reference`. Only rows with `TRUE` can enter the active effect
table or Layer 1 penalty projection.

## Stage 3 versus Stage 4 use

Stage 3 uses active condition coefficients `beta[e, c]` to identify supported
metabolic target genes and complete-GPR core reactions. It does not require a
reference contrast.

Stage 4 uses only comparable condition effects `Delta beta[e, c]`. For metacell
`u`, RegCompass reconstructs the Pando model space and computes:

\[
M_{g,c,u}=\sqrt{clamp(R^2_{g,c},0,1)}
\sum_{e:target(e)=g}\Delta\beta_{e,c}z_{e,u}.
\]

The signed bounded modifier is:

\[
R_{g,c,u}=\tanh(M_{g,c,u}).
\]

No metacell-wise robust rescaling or absolute-sum coefficient normalization is
applied. Doubling all comparable coefficient effects therefore doubles the raw
projection.

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
pando_objects/condition_grn_fit_v2.rds
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
