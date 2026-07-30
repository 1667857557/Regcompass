# Stage input-output contracts

RegCompass connects stages only when classes, workflow settings, GEM provenance,
analysis mode, metacell construction provenance, and unit order agree.

## Stage 1: Pando evidence

Class: `regcompass_grn_step`.

Common outputs:

```r
step1$grn_result$analysis_mode
step1$grn_result$condition_coefficients_calculated
step1$grn_result$target_metabolic_genes
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$params$requested_condition_col
step1$params$condition_col
step1$params$celltype_col
```

### `analysis_mode = "condition_grn"`

Selected only when the effective condition column contains at least two levels.
Pando uses the canonical unversioned `pando_condition_grn_fit` schema, aligned
absolute coefficients, equal-condition transforms, nested outer-heldout
projections, estimability/support masks, and exactly-once OOF assignment.

Additional output:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition_effect
```

### `analysis_mode = "standard_pando"`

Selected when the condition column is omitted, absent, or contains one level.
RegCompass calls original `Pando::infer_grn()` independently within each broad
cell type. No `ConditionGRNFit`, condition coefficient, condition deviation, or
condition contrast is calculated.

```r
step1$grn_result$standard_pando_objects
step1$grn_result$condition_grn_fits       # empty list
step1$grn_result$condition_coefficients_calculated  # FALSE
```

The effective constant condition label is retained only so downstream tables
have one grouping value. It is not used to fit the standard Pando model.

## Stage 2: native SuperCell metacells

Class: `regcompass_metacell_step`.

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$metacell_object
step2$params
```

RegCompass passes the two biological labels directly to SuperCell:

```r
SuperCell::SCimplify_from_embedding(
  X = joint_embedding,
  cell.annotation = cell_type,
  cell.split.condition = condition,
  gamma = gamma
)
```

For standard mode with an omitted condition,
`cell.split.condition = NULL`. RegCompass does not create or pass a concatenated
`condition__cell_type` stratum. The cache schema is
`regcompass_native_supercell_metacell_cache_v1` and records:

```text
native_supercell_api = SCimplify_from_embedding
condition_argument = cell.split.condition
celltype_argument = cell.annotation
temporary_combined_stratum = FALSE
```

RNA and ATAC raw counts are aggregated from the exact returned
`membership(cell_id, metacell_id)` table.

## Stage 3: biological meta-modules

Class: `regcompass_meta_module_step`.

Active standard or condition-aware Pando target genes form one supported gene
set per effective condition and broad cell type. Positive and negative
coefficients both count as regulatory evidence. A reaction is a core only when
one complete GPR branch is contained in the supported set.

Expansion is one ordered pass:

1. core subsystem;
2. direct KEGG/Reactome equivalence;
3. direct master-Rhea equivalence.

Stage 3 creates a reaction catalogue, not a GEM, and does not run FASTCORE.

## Stage 4: regulatory Layer 1

Class: `regcompass_layer1_step`.

Both modes use cell-first TF RNA × peak ATAC projections followed by exact
SuperCell membership aggregation. Interactions are never reconstructed from
metacell means.

- condition mode: outer-heldout common-support Pando projections;
- standard mode: original Pando full-fit coefficients, with no condition
  coefficient calculation.

The modes then share the same processing:

```text
cell-level regulatory projection
→ metacell mean
→ cell-type Gamma–Poisson latent RNA support
→ reliability × tanh(projection / shared scale)
→ bounded RNA-support odds modifier
→ GPR reaction expression
```

`regulatory_alpha` is fixed at `1`. A non-finite target modifier is neutralized
to `R = 0`, giving exactly RNA-only support for that gene–metacell entry.

## Stage 5: Layer 2

Class: `regcompass_layer2_step`.

For each medium, one shared union GEM and one global FASTCORE completion are used
for all metacells and both evidence routes. GPR OR branches are summed while
unavailable branches are ignored. Missing final reaction expression is assigned
`E = 0` before conversion and therefore receives expression-linked penalty `1`.

Required outputs include:

```r
step5$penalty
step5$vmax
step5$score
step5$model_cache_summary
step5$structural_model_contract
```

## Stage 6: results

The final result records the selected mode explicitly:

```r
result$analysis_mode
result$condition_coefficients_calculated
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

For one effective condition, `reaction_ranking` and `condition_summary` are
returned and `condition_contrast` is empty. No artificial second condition or
condition coefficient is generated.
