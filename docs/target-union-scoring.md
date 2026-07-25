# Direct database-linked non-core scoring in existing union GEMs

`rc_regcompass_step_target_union()` is an optional second-pass analysis after a completed stepwise Layer 2 run with `model_mode = "meta_module_gem"`.

## Structural source

The function reuses the exact medium-specific union GEM files recorded by:

```r
layer2$model_cache_summary
```

These models already contain:

- the merged biological meta-module reactions;
- medium-specific bounds;
- global FASTCORE support;
- the original direction-specific core target set.

The Stage 3 object supplies anchor provenance through:

```r
meta_modules$merged_modules$merged_core_reactions
meta_modules$merged_modules$merged_reaction_membership
```

These Stage 3 tables are not GEMs.

## Mapping rule

Selected previous core reactions are used as anchors. Direct non-core targets are reactions sharing at least one:

- KEGG reaction identifier;
- Reactome reaction identifier;
- master Rhea identifier.

No subsystem expansion, transitive expansion, metabolite-neighbour expansion, or new FASTCORE completion is performed.

## Availability rule

A mapped reaction is scoreable only when it is present in the actual cached union GEMs required for the restart. This allows globally added FASTCORE support reactions to be scored when they are present in all reused models, even if they were absent from the Stage 3 merged biological catalogue.

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
targeted$previous_union_membership
targeted$microcompass
```

The output file containing reused model membership is:

```text
reused_union_gem_membership.tsv.gz
```

Removed Stage 3 names such as `global_modules`, `global_core_reactions`, and `global_reaction_membership` must not be used.
