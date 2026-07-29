# Tutorial 4: targeted reaction scoring

Use `rc_regcompass_step_target_union()` after a completed
`model_mode = "meta_module_gem"` run to score direct database equivalents of
selected reaction anchors.

This function reuses the cached Stage 5 model. It does not rebuild the model or
rerun FASTCORE.

## Load stages

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

`step5$model_cache_summary$file` must point to an available completed model.

## Select reaction anchors

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted/reaction_anchors",
  core_reaction_ids = c("MAR04381", "MAR04379", "MAR04391"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

Despite the retained argument name, `core_reaction_ids` may contain an original
core reaction or another reaction present in the supplied GEM. The anchor is
used to find direct KEGG, Reactome, or master-Rhea equivalents and is not
automatically rescored.

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
    solver = "highs"
  )
)
```

- `complete_gpr`: selected genes must satisfy a complete GPR branch in the original core set;
- `any_direct`: match any directly associated gene in the original core set.

## Mapping scope

The function includes only direct equivalents sharing a KEGG reaction ID,
Reactome reaction ID, or master Rhea ID with an anchor.

It does not perform subsystem, transitive, metabolite-neighbour, or one-hop
expansion. A mapped target is scored only when it is available in every required
cached model.

## Inspect outputs

```r
targeted$selected_anchor_reactions
targeted$selected_core_reactions
targeted$selected_noncore_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$microcompass$model_cache_summary
```

Relation-level provenance:

```r
targeted$expanded_reaction_catalog[, c(
  "anchor_reaction_id",
  "anchor_is_original_core",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "available_in_all_cached_union_gems"
)]
```

Target-level output:

```r
targeted$expanded_scoring_targets[, c(
  "reaction_id",
  "anchor_reaction_ids",
  "expansion_types",
  "source_annotations"
)]
```

The persistent catalogue is written to
`merged_meta_module_catalogue_membership.tsv.gz`. The model cache checksum is
retained to verify exact structural reuse.

Public API: [functions.md](functions.md). Mathematical definitions:
[Mathematical model](mathematical-model.md).
