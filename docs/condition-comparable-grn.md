# Shared-design, condition-comparable Pando GRNs

RegCompass Stage 1 passes the complete normalized multiome object once through
`Pando::initiate_grn()` and `Pando::find_motifs()`, then calls
`Pando::infer_condition_grn()`. Pando fits one cell-type design across its
conditions and returns a versioned `ConditionGRNFit`.

The current public entry point is `rc_regcompass_step_grn()`. The complete-run
wrappers call the same function and preserve the same fit contract.

## Why the default is shared-design independent

The objective is to approximate separate condition-specific Pando fits while
making their coefficients directly comparable. The default
`shared_design_independent` engine therefore shares:

- the complete TF–peak–target edge dictionary;
- the edge-by-condition eligibility mask;
- the pooled transform of every final `TF RNA × peak ATAC` predictor;
- the target transform;
- within each target, the lambda path and selected lambda.

It does not add a group penalty that shrinks the same edge across conditions.
For condition \(c\), Pando minimizes an ordinary elastic-net objective at the
target-specific \(\lambda\) shared by conditions:

\[
\frac{1}{2n_c}\lVert y_c-\alpha_c-X_c\beta_c\rVert_2^2
+\lambda\left[
  \eta\lVert\beta_c\rVert_1+
  \frac{1-\eta}{2}\lVert\beta_c\rVert_2^2
\right].
\]

The condition coefficient vectors are separable at fixed \(\lambda\).
Cross-condition information is used only to define the common design,
transform, validation loss, and tuning choice.

RegCompass enforces:

```r
pando_infer_args = list(
  method = "shared_design_independent",
  candidate_screen = "condition_union",
  condition_mix = 1,
  condition_weight = "equal",
  scale = TRUE
)
```

`reference_condition` should be set explicitly for a stable biological
interpretation, for example `"Control"`. When omitted, Pando uses the first
condition level within each cell type.

## Exact edge union

For `candidate_screen = "condition_union"`, the TF and peak of an edge must
both pass screening inside the same condition. Pando then takes the union of
those complete edges:

\[
\mathcal E =
\bigcup_c
\{(TF,peak,target): TF\text{ and }peak\text{ pass in }c\}.
\]

This prevents a false edge from being created by pairing a TF retained only in
condition A with a peak retained only in condition B. The
`eligibility_mask[edge, condition]` records the exact result. A structurally
ineligible coefficient is fixed at zero and remains distinguishable from an
eligible fitted zero.

## Common coefficient units

Pando first forms the raw edge predictor

\[
x_{e,u} =
\mathrm{RNA}_{TF(e),u}\,
\mathrm{ATAC}_{peak(e),u}.
\]

It then calculates one center \(\mu_e\) and scale \(s_e\) across all cells of
the cell type and uses

\[
z_{e,u} = \frac{x_{e,u}-\mu_e}{s_e}.
\]

Scaling the final interaction is important: separately scaling RNA and ATAC
before multiplication does not produce the same predictor. The target is also
centered and scaled once over the pooled cell type.

The complete values are retained in:

- `edge_table`
- `beta`
- `contrast`
- `eligibility_mask`
- `predictor_transform`
- `response_transform`
- `target_fit`
- `target_rsq`

RegCompass writes the fit contracts to
`pando_condition_grn_fits.rds` and the edge transforms to
`pando_edge_predictor_transforms.tsv.gz`.

The in-memory current API is:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
```

Historical `sample_status`, `tf_peak_gene_all`, and
`tf_peak_gene_significant` fields are retained only as compatibility aliases.

## Reference contrasts

Condition effects use an explicit reference:

\[
\Delta\beta_{e,c}
=\beta_{e,c}-\beta_{e,reference}.
\]

The reference column is exactly zero. The Universal Pando `Network` remains an
equal-condition coefficient mean for compatibility and visualization only; it
is never used as the contrast baseline.

Stage 3 uses active \(\beta_{e,c}\) values to identify supported target genes
and complete-GPR core reactions. The condition-effect table
`tf_peak_gene_condition_effect` contains \(\Delta\beta\) and is used only for
the downstream regulatory modifier.

## Exact RegCompass model-space projection

For metacell \(u\), RegCompass reads the same TF RNA and peak ATAC features,
applies Pando's stored transform, and computes

\[
M_{g,c,u} =
\sqrt{\operatorname{clamp}(R^2_{g,c},0,1)}
\sum_{e:\,target(e)=g}
\Delta\beta_{e,c}\,z_{e,u}.
\]

No metacell-wise robust rescaling is applied. Edge weights are not divided by
\(\sum_e|\Delta\beta_e|\); doubling every fitted effect therefore doubles the
raw projection instead of disappearing under normalization.

The bounded signed modifier is

\[
R_{g,c,u}=\tanh(M_{g,c,u}).
\]

It updates bounded target-gene RNA support \(C\) on the log-odds scale:

\[
C_{\mathrm{multiome}} =
\frac{C\,2^{\alpha R}}
{1-C+C\,2^{\alpha R}}.
\]

Layer 1 stores both:

- `gene_regulatory_model_projection`: raw \(M\);
- `gene_regulatory_modifier`: bounded \(R\).

## One shared metabolic model

After regulatory evidence defines the merged biological reaction catalogue,
RegCompass applies each medium and performs one global FASTCORE completion.
All conditions and metacells in that medium reuse the same stoichiometric
matrix and identical lower/upper bounds. Condition differences enter the
directional scoring problem only through the reaction penalty derived from
Layer 1 evidence.

## Output audit

Stage 1 writes:

- `pando_group_status.tsv.gz`
- `pando_tf_peak_gene_condition_all.tsv.gz`
- `pando_tf_peak_gene_condition_active.tsv.gz`
- `pando_tf_peak_gene_condition_effect_all.tsv.gz`
- `pando_tf_peak_gene_condition_effect_active.tsv.gz`
- `pando_tf_peak_gene_universal.tsv.gz` (summary only)
- `pando_condition_network_index.tsv.gz`
- `pando_condition_fit_diagnostics.tsv.gz`
- `pando_edge_predictor_transforms.tsv.gz`
- `pando_condition_grn_fits.rds`
- `pando_objects/condition_grn_fit_v2.rds`, when enabled

These artifacts preserve the inputs required to reproduce the coefficient
contrast and RegCompass projection without re-fitting a GRN.
