# Tutorial Level 4: remap selected genes or reactions and run second-pass scoring

Use this tutorial after a completed stepwise `meta_module_gem` analysis when selected complete-GPR core reactions should be used as anchors to identify and score directly database-linked non-core reactions.

This operation reuses the **exact final medium-specific union GEM files** created by Stage 5. It does not treat the Stage 3 merged meta-module catalogue as a GEM.

## Load the completed stages

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

Required conditions:

```r
stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  inherits(step4, "regcompass_layer1_step"),
  inherits(step5, "regcompass_layer2_step"),
  identical(step5$model_mode, "meta_module_gem"),
  nrow(step5$model_cache_summary) > 0,
  all(c(
    "medium_scenario",
    "file",
    "file_checksum",
    "build_strategy",
    "completion_stage"
  ) %in% colnames(step5$model_cache_summary))
)
```

## Select anchors by reaction ID

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted/reaction_anchors",
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

Selected reaction IDs must have been complete-GPR core targets in the original Stage 5 run.

## Select anchors by gene

```r
targeted_gene <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted/gene_anchors",
  core_genes = c("SLC7A11", "GCLC", "GCLM"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 600
  )
)
```

`gene_match = "complete_gpr"` requires the selected genes to satisfy at least one full GPR isozyme group. Use `"any_direct"` only for exploratory mapping.

## Mapping scope

The second-pass target catalogue includes non-core reactions sharing a direct:

- KEGG reaction identifier;
- Reactome reaction identifier;
- master Rhea identifier

with one or more selected core anchors.

It does not perform subsystem expansion, transitive cross-reference expansion, metabolite-neighbour expansion, model reconstruction, or FASTCORE completion.

## Exact model-reuse contract

A mapped reaction is scoreable only when it is present in every final medium-specific union GEM required by the analysis. For each Stage 5 cache row, RegCompass validates:

- the cached file path;
- `file_checksum`;
- `build_strategy = "medium_specific_union_gem"`;
- `completion_stage = "single_global_fastcore_after_meta_module_merge"`;
- `model$is_union_gem`;
- the model's medium-scenario identifier.

The exact stoichiometric matrix, bounds, medium constraints, and support reactions are reused. The second pass does not rebuild a GEM and does not rerun FASTCORE.

```r
targeted$selected_core_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$merged_catalogue_membership
targeted$microcompass$model_cache_summary

targeted$microcompass$params[c(
  "structural_model_reused_exactly",
  "fastcore_rerun",
  "model_rebuild"
)]
```

The expected flags are:

```r
stopifnot(
  isTRUE(targeted$microcompass$params$structural_model_reused_exactly),
  identical(targeted$microcompass$params$fastcore_rerun, FALSE),
  identical(targeted$microcompass$params$model_rebuild, FALSE)
)
```

## Inspect provenance

Relation-level provenance:

```r
targeted$expanded_reaction_catalog[, c(
  "anchor_core_reaction_id",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "present_in_merged_catalogue",
  "merged_catalogue_inclusion_stage",
  "available_in_all_cached_union_gems"
)]
```

Reaction-level aggregated targets:

```r
targeted$expanded_scoring_targets[, c(
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations",
  "merged_catalogue_inclusion_stage"
)]
```

`merged_catalogue_inclusion_stage` describes whether the reaction was already a biological member of the merged Stage 3 catalogue. Reactions absent from that catalogue can still be scoreable when global FASTCORE added them to every reused final union GEM.

The persistent catalogue file is:

```text
merged_meta_module_catalogue_membership.tsv.gz
```
