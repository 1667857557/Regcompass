# Tutorial Level 4: remap selected genes or reactions and run second-pass scoring

Use this tutorial after a completed stepwise `meta_module_gem` analysis when selected complete-GPR core reactions should be used as anchors to identify and score directly database-linked non-core reactions.

This operation reuses the **actual medium-specific union GEM files** created by Stage 5. It does not treat the Stage 3 merged meta-module catalogue as a GEM.

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
  nrow(step5$model_cache_summary) > 0
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

It does not perform subsystem expansion, transitive cross-reference expansion, metabolite-neighbour expansion, or a new FASTCORE reconstruction.

## Why no FASTCORE is rerun

The target reaction must already be present in every reused medium-specific union GEM required by the analysis. The exact cached model structure and medium bounds are reused. This keeps the second-pass score directly comparable with the original Layer 2 structural context.

```r
targeted$selected_core_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$microcompass$model_cache_summary
```

## Inspect provenance

```r
targeted$expanded_reaction_catalog[, c(
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations",
  "previous_union_inclusion_stage"
)]
```

`previous_union_inclusion_stage` describes the reaction's status in the Stage 5 union GEM context. The Stage 3 catalogue remains available separately as:

```r
step3$merged_modules$merged_reaction_membership
```

Do not use the removed fields `global_modules`, `global_core_reactions`, or `global_reaction_membership`.
