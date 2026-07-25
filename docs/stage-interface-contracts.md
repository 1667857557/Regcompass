# Stage input-output contracts

RegCompassR 1.8.4 connects stages only when classes, workflow settings, GEM provenance, metacell construction provenance, and scoring-unit order agree.

## Stage 1: Pando evidence

Class: `regcompass_grn_step`

Required outputs:

```r
step1$grn_result$target_metabolic_genes
step1$grn_result$tf_peak_gene_significant
step1$grn_result$sample_status
step1$grn_result$normalization_policy$pando_regions
step1$gem_fingerprint
step1$params
```

Every scored `condition × cell type` group must have a successful Pando fit with significant edges. `target_metabolic_genes` is the intersection of GEM GPR genes and RNA-assay row names. For human analyses without an explicit `pando_initiate_args$regions`, `pando_regions` records the union of Pando's hg38 phastCons and SCREEN ccRE data objects.

## Stage 2: metacells

Class: `regcompass_metacell_step`

Required outputs:

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$metacell_object
step2$params
```

The merged metacell object and metadata must contain the same ordered units. Reduction names, dimensions, cell labels, assay fingerprints, and embedding fingerprints are part of the cache contract.

## Stage 3: biological meta-modules

Class: `regcompass_meta_module_step`

Required outputs:

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$group_coverage
```

Contract:

- `supported_metabolic_genes` contains one row per condition, cell type, and Human-GEM target gene with at least one significant Pando TF–peak–gene row;
- positive and negative Pando coefficients both count as regulatory evidence;
- all supported genes in one `condition × cell type` form one GPR-evaluation set;
- Stage 3 performs no shared-TF target projection, top-k graph pruning, or connected-component analysis;
- `core_gene_reaction` marks a reaction as core only when one complete GPR branch is contained in the supported gene set;
- `merged_core_reactions` contains deduplicated complete-GPR core reactions;
- `merged_reaction_membership` contains deduplicated biological reactions only;
- `merged_modules$is_gem` is `FALSE`;
- `merged_modules$fastcore_applied` is `FALSE`;
- Stage 3 does not apply medium constraints or run FASTCORE;
- the merged object is a catalogue and must not be described as a union GEM.

The only Stage 3 parameters are:

```r
meta_module_args = list(
  subsystem_table = NULL,
  expansion_mode = "ordered_once",
  max_iterations = 10
)
```

## Stage 4: Layer 1

Class: `regcompass_layer1_step`

Required outputs:

```r
step4$reaction_expression
step4$metacell_meta
step4$workflow_params
step4$gem_fingerprint
```

The reaction-expression matrix must contain every merged core reaction and the same ordered metacells represented by Stage 2.

## Stage 5: Layer 2

Class: `regcompass_layer2_step`

For `model_mode = "meta_module_gem"`, required inputs are:

```r
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

Required outputs include:

```r
step5$model_cache_summary
step5$source_core_reactions
step5$source_merged_reaction_membership
step5$score
step5$penalty
step5$vmax
```

Each `model_cache_summary` row identifies one final medium-specific union GEM and records at least:

```r
medium_scenario
file
file_checksum
build_strategy
completion_stage
```

The cached model must record:

```r
model$is_union_gem
model$union_gem_medium_scenario
model$build_params$completion_stage
model$reaction_meta$global_fastcore_support
```

The required build metadata are:

```text
build_strategy = medium_specific_union_gem
completion_stage = single_global_fastcore_after_meta_module_merge
```

All conditions and metacells within the same medium must resolve to the same model file. Different media may resolve to different union-GEM structures because their global FASTCORE support sets may differ.

## Stage 6: results

The final result contains:

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$merged_grn_meta_modules
result$microcompass
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

`merged_grn_meta_modules` is the Stage 3 catalogue. `microcompass$model_cache_summary` identifies the final Stage 5 union GEMs.

## Global FASTCORE configuration

The only structural completion controls are supplied at Stage 5:

```r
layer2_args = list(
  model_params = list(
    completion_time_limit = 600,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE
  )
)
```

## Target-union restart contract

`rc_regcompass_step_target_union()` requires:

- the original Stage 3 merged catalogue for anchor provenance;
- the original Layer 1 matrix;
- the completed `meta_module_gem` Stage 5 object;
- accessible final Stage 5 union-GEM files;
- matching model-file checksums and medium identifiers.

The selected genes or reaction IDs determine mapping anchors only. Target availability and all LP calculations are evaluated in the exact cached final union GEMs. The second pass does not rebuild a GEM, does not change medium bounds, and does not rerun FASTCORE.
