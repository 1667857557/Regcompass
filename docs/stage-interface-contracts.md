# Stage input-output contracts

RegCompassR 1.8.4 connects stages only when classes, workflow settings, GEM provenance, metacell construction provenance, and scoring-unit order agree.

## Stage 1: GRN

Class: `regcompass_grn_step`

Required outputs:

```r
step1$grn_result$tf_peak_gene_significant
step1$grn_result$sample_status
step1$gem_fingerprint
step1$params
```

Every scored `condition × cell type` group must have a successful GRN with significant edges.

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
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$group_coverage
```

Contract:

- `merged_core_reactions` contains deduplicated complete-GPR core reactions;
- `merged_reaction_membership` contains deduplicated biological reactions only;
- `merged_modules$is_gem` is `FALSE`;
- `merged_modules$fastcore_applied` is `FALSE`;
- Stage 3 does not apply medium constraints or run FASTCORE.

Removed fields:

```text
global_modules
global_core_reactions
global_reaction_membership
local_completed_reaction_membership
local_fastcore_summary
local_fastcore_diagnostics
local_fastcore_completion_iterations
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

Each `model_cache_summary` row identifies one cached medium-specific union GEM. Its model must record:

```r
model$is_union_gem
model$union_gem_medium_scenario
model$build_params$completion_stage
model$reaction_meta$fastcore_support
```

The only permitted completion stage is:

```text
single_global_fastcore_after_meta_module_merge
```

All conditions and metacells within the same medium must resolve to the same model file.

## Stage 6: results

The final result contains:

```r
result$condition_grn_meta_modules
result$merged_grn_meta_modules
result$microcompass
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

`global_grn_meta_modules` is removed. `grn_meta_modules` remains a generic alias of `merged_grn_meta_modules`.

## Canonical runner argument contract

The following are rejected:

```r
layer1_args = list(local_fastcore = TRUE)
layer1_args = list(local_fastcore_args = list(...))
```

Use:

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

- the original Stage 3 merged catalogue;
- the original Layer 1 matrix;
- the completed `meta_module_gem` Stage 5 object;
- accessible cached union-GEM files.

Target availability is checked against the actual cached union GEMs, not against the pre-FASTCORE Stage 3 catalogue alone.
