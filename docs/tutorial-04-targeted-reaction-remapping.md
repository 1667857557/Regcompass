# Tutorial 4: targeted reaction remapping

Use `rc_regcompass_step_target_union()` after a completed
`model_mode = "meta_module_gem"` run to score direct database equivalents of
selected reaction anchors.

This is an optional second LP pass, not a condition-comparability guardrail. It
reuses the exact cached Stage 5 model, does not rebuild the GEM, and does not
rerun FASTCORE. In condition-aware runs, the Layer 1 `reaction_expression`
input is the primary BH-filtered fixed-dictionary regulatory route. Historical
fields containing `condition_full_oof` remain compatibility aliases and do not
indicate OOF estimation.

Mathematical definitions are in
[Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

## Load completed stages

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

`step5$model_cache_summary$file` must point to an available completed model. The
function validates the stored checksum and union-GEM provenance before scoring.

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
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

- `complete_gpr`: selected genes must satisfy a complete GPR branch in the
  original core set;
- `any_direct`: match any directly associated gene in the original core set.

## Mapping scope

The function includes only direct equivalents sharing a KEGG reaction ID,
Reactome reaction ID, or master Rhea ID with an anchor.

It does not perform subsystem, transitive, metabolite-neighbour, or one-hop
expansion. A mapped target is scored only when it is available in every required
cached medium-specific union GEM. Reactions already scored as original Layer 2
cores are not recomputed.

## Regulatory evidence used by the second pass

The second pass uses:

```r
step4$reaction_expression
```

For a condition-aware run this is the primary fixed-dictionary condition route.
The historical field below is an equal-valued compatibility alias:

```r
step4$reaction_expression_condition_full_oof
```

The targeted reactions therefore use the same condition-specific regulatory
evidence as the original Stage 5 ranking. Targets without significant estimable
condition-GRN edges retain the same RNA-only neutral fallback used in Stage 4.

## Structural reuse

The following are fixed and reused:

- medium-specific union GEM file and checksum;
- reaction order and bounds;
- global FASTCORE completion;
- target direction and `omega`;
- medium-specific `vmax` definition;
- metacell order and condition/cell-type metadata;
- Layer 1 GPR aggregation settings.

The targeted pass does not alter the original Stage 5 cache.

## Parallel execution

The same Layer 2 backend can be reused:

```r
library(BiocParallel)

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 30L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 30L, progressbar = TRUE)
}
```

Each worker forces numerical libraries and HiGHS to one internal thread. Reduce
the number of workers when the cached GEM or target set makes per-worker memory
the limiting resource.

## Inspect outputs

```r
targeted$selected_core_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
targeted$merged_catalogue_membership
targeted$microcompass$model_cache_summary
targeted$microcompass$penalty
targeted$microcompass$score
```

Relation-level provenance:

```r
targeted$expanded_reaction_catalog[, c(
  "anchor_core_reaction_id",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "present_in_merged_catalogue",
  "available_in_all_cached_union_gems"
)]
```

Target-level output:

```r
targeted$expanded_scoring_targets[, c(
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations"
)]
```

The persistent catalogue is written to
`merged_meta_module_catalogue_membership.tsv.gz`. The model cache checksum is
retained to verify exact structural reuse.

## When to rerun earlier stages

Changing only the selected anchors does not require rerunning Stages 1–5.
Changing the Stage 4 reaction evidence, Stage 5 medium, GEM, FASTCORE completion,
bounds, direction, or target-flux controls invalidates the targeted pass and
requires regenerating the relevant upstream cache first.


## Cell-type-scoped cache reuse

Targeted remapping preserves the original structural boundary. A target linked
to a core reaction from cell type `c` is evaluated only in cached
`cell_type = c` union GEMs. The function intersects candidate reactions across
media within that cell type; it never intersects or merges cached GEMs from
different cell types and never reruns FASTCORE.
