# Pando common-dictionary condition-GRN contract

This page defines the implemented Pando–RegCompass interface for paired
single-cell RNA+ATAC data.

## Automatic routing

```text
at least two condition levels → common-dictionary condition GRNs
zero or one condition level   → original per-cell-type Pando GRN
```

`rc_regcompass_step_grn()` always fits broad cell types independently. Sample or
donor metadata are retained only as provenance and are not fitting gates.

## Multiple-condition algorithm

For one broad cell type and target gene `g`, candidate discovery is performed on
the complete eligible cell type and separately in every condition:

\[
E_g^{global}=Discover(D_{all},g),
\]

\[
E_g^{(c)}=Discover(D_c,g).
\]

`Discover()` contains the original Pando biological candidate steps:

- peak-to-gene regulatory domain;
- motif-to-TF mapping;
- peak-target correlation threshold;
- TF-target correlation threshold.

The final dictionary is the exact triple union:

\[
E_g^\cup=E_g^{global}\cup\bigcup_c E_g^{(c)}.
\]

The union key is `(target, TF, region)`. TF, peak and target node sets are not
unioned separately and are never recombined by Cartesian product.

After the union is frozen, each condition fits exactly the same predictor set:

\[
Y_{igc}=\alpha_{gc}+\sum_{r\in E_g^\cup}\beta_{rgc}
E_{i,t(r),c}A_{i,p(r),c}+\varepsilon_{igc}.
\]

The canonical controls are:

```r
pando_infer_args = list(
  tf_cor = 0.1,
  peak_cor = 0,
  adjust_method = "BH",
  padj_threshold = 0.05,
  rank_action = "mark",
  min_residual_df = 1L
)
```

The final model is Gaussian identity GLM with `interaction = ":"` and
`scale = FALSE`. RNA and ATAC must therefore be preprocessed before condition
splitting. All conditions use the same assay layers, variable definitions,
coefficient units and conditioning covariates. No global coefficient is used to
rescale or calibrate a condition coefficient.

## Effect and significance contract

For edge `r` and condition `c`:

\[
condition\_effect_{r,c}=\widehat\beta_{r,c}.
\]

Pando stores:

```text
estimate
std_err
statistic
pval
padj
significant
penalty_effect
direction
estimable
zero_variance
aliased
```

`estimate` is always the complete fitted coefficient when estimable.
`penalty_effect` is defined as:

\[
penalty\_effect_{r,c}=
\begin{cases}
\widehat\beta_{r,c}, & padj_{r,c}<0.05,\\
0, & otherwise.
\end{cases}
\]

edge enters the active penalty table. These gates do not rewrite the complete
coefficient table.

A zero-variance or aliased edge has `estimate = NA` and `estimable = FALSE`. It
is not interpreted as an estimated biological zero. Its realized downstream
penalty contribution is neutral because no finite `penalty_effect` exists.

The ordinary GLM P values are conditional on the frozen candidate dictionary.
They do not include selective-inference correction for the preceding correlation
screening.

## Why coefficients are comparable

Comparability follows from four simultaneous constraints:

```text
one global preprocessing reference
+ one exact target-specific edge dictionary
+ one Gaussian identity interaction formula
+ scale = FALSE
```

Thus each condition coefficient has the same variable definition and numerical
unit. Conditions may retain different significant edge sets and may estimate
opposite signs on the same edge, but the coefficients are evaluated on the same
partial-regression coordinate.

## Pando fit schema

The canonical schema is:

```text
pando_condition_grn_common_dictionary_v1
```

Each `ConditionGRNFit` contains:

```text
cell_type
condition_levels
condition_cell_ids
edge_dictionary
coefficients
fit
network_names
padj_threshold
adjust_method
scale
interaction
projection_effect_column
projection_policy
```

Every condition also has a standard Pando `Network` object, so `coef()`, `gof()`
and existing Pando network utilities remain available.

## RegCompass Stage 1 artifacts

```text
pando_group_status.tsv.gz
pando_tf_peak_gene_condition_all.tsv.gz
pando_tf_peak_gene_condition_active.tsv.gz
pando_tf_peak_gene_universal.tsv.gz
pando_condition_grn_fits.rds
```

`condition_all` is the lossless coefficient table. `condition_active` contains
only effects eligible for penalty after significance, effect-size and target
model-fit gates.

## Paired-cell projection and metacell aggregation

For a paired cell `i` and target `g`, Pando calculates:

\[
G_{igc}=\sum_{r\in E_g^\cup}
penalty\_effect_{r,c}E_{i,t(r),c}A_{i,p(r),c}.
\]

RegCompass then averages `G` within the exact SuperCell membership. Projection
must occur before aggregation. The workflow does not reconstruct TF×ATAC from
metacell averages, refit coefficients after aggregation, or normalize each
condition separately.

## Stage 4 and Stage 5 compatibility fields

The existing public result schema retains historical names containing `_oof`,

```text
gene_projection = primary BH-filtered full-fit projection
penalty          = primary metabolic penalty
```

These names no longer denote OOF estimation or a jointly estimable common-support
model.

## Retired parameters and runtime

The following condition-GRN controls are rejected:

```text
candidate_screen
condition_mix
condition_weight
alpha
nlambda / lambda / lambda_min_ratio
outer_nfolds / inner_nfolds
lambda_selection
scale
engine_control
comparison_conditions
```

runtime is not part of the current Pando package.

## No-condition and single-condition behavior

If `condition_col` is not supplied, is absent, or has fewer than two observed
levels, RegCompass fits the original Pando Gaussian interaction GRN independently
for each retained broad cell type. No `ConditionGRNFit`, condition coefficient,
condition contrast or synthetic condition projection is created.
