# Tutorial Level 4: remap selected genes or reactions

Use this tutorial after a completed stepwise `meta_module_gem` analysis to score non-core reactions that are directly linked to selected complete-GPR core reactions.

The second pass reuses the cached Stage 5 models. It validates the cache checksums and medium identities, and does not rebuild the model or rerun FASTCORE.

## Load the completed stages

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

`step5` must come from a completed `model_mode = "meta_module_gem"` run with available model-cache files.

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

The second pass includes non-core reactions that share a direct KEGG reaction ID, Reactome reaction ID, or master Rhea ID with a selected core anchor.

It does not perform subsystem, transitive, metabolite-neighbour, or one-hop expansion. A mapped reaction is scored only when it is present in every required cached Stage 5 model.

## Inspect outputs and provenance

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

Reaction-level targets:

```r
targeted$expanded_scoring_targets[, c(
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations",
  "merged_catalogue_inclusion_stage"
)]
```

The persistent mapping table is written to:

```text
merged_meta_module_catalogue_membership.tsv.gz
```
