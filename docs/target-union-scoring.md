# Direct database-linked non-core scoring in final union GEMs

`rc_regcompass_step_target_union()` is an optional second-pass analysis after a completed stepwise Layer 2 run with `model_mode = "meta_module_gem"`.

## Structural source

The function reuses the exact final medium-specific union GEM files recorded by:

```r
layer2$model_cache_summary
```

Each cache row must record the model file, checksum, medium scenario, build strategy, and completion stage. The cached models already contain:

- the merged biological meta-module reactions;
- medium-specific bounds;
- global FASTCORE support;
- the original direction-specific core target set.

The Stage 3 object supplies anchor provenance through:

```r
meta_modules$merged_modules$merged_core_reactions
meta_modules$merged_modules$merged_reaction_membership
```

These Stage 3 tables are catalogue tables, not GEMs.

## Mapping rule

Selected merged core reactions are used as anchors. Direct non-core targets are reactions sharing at least one:

- KEGG reaction identifier;
- Reactome reaction identifier;
- master Rhea identifier.

No subsystem expansion, transitive expansion, metabolite-neighbour expansion, model reconstruction, or FASTCORE completion is performed.

## Availability and validation rule

A mapped reaction is scoreable only when it is present in every required final union GEM. Before scoring, RegCompass verifies:

```text
file_checksum
build_strategy = medium_specific_union_gem
completion_stage = single_global_fastcore_after_meta_module_merge
model$is_union_gem = TRUE
model$union_gem_medium_scenario matches the cache row
```

This allows globally added FASTCORE support reactions to be scored when they are present in all reused final models, even when absent from the Stage 3 merged biological catalogue.

## Example

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted",
  core_reaction_ids = c("MAR04324"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 600
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

## Outputs

```r
targeted$selected_core_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$merged_catalogue_membership
targeted$microcompass
```

The scoring result records the exact reuse policy:

```r
targeted$microcompass$params[c(
  "structural_model_reused_exactly",
  "fastcore_rerun",
  "model_rebuild"
)]
```

The output catalogue file is:

```text
merged_meta_module_catalogue_membership.tsv.gz
```

The expanded tables use:

```text
present_in_merged_catalogue
merged_catalogue_is_core
merged_catalogue_inclusion_stage
available_in_all_cached_union_gems
```
