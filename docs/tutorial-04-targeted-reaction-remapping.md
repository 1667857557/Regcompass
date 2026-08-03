# Tutorial 4: targeted reaction remapping

Use `rc_regcompass_step_target_union()` after a completed
`model_mode = "meta_module_gem"` run to score direct database equivalents of
selected reaction anchors.

This is an optional second LP pass, not a condition-comparability guardrail. It
reuses the exact cached Stage 5 models for the corresponding cell types and
media, does not rebuild a GEM, and does not rerun FASTCORE. In condition-aware
runs, the Layer 1 `reaction_expression` input is the primary BH-filtered
fixed-dictionary regulatory route. Historical fields containing
estimation.

Mathematical definitions are in
[Tutorial 3](tutorial-03-mathematical-model.md). Public API:
[functions.md](functions.md).

## Load completed stages

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

`step5$model_cache_summary` must contain an available model for each retained
`cell_type × medium_scenario` combination. The targeted function validates the
stored file, checksum, cell type, medium, build strategy, FASTCORE completion
stage, and union-GEM provenance before scoring.

Inspect the structural source first:

```r
step5$model_cache_summary[, c(
  "cell_type",
  "medium_scenario",
  "file",
  "file_checksum",
  "n_celltype_biological_reactions",
  "n_celltype_fastcore_support_reactions"
)]
```

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
core reaction or another reaction present in the supplied GEM. Each anchor is
evaluated separately in the cell types where it belongs to the original core or
is available in the corresponding cached structural model. The anchor is used
to find direct KEGG, Reactome, or master-Rhea equivalents and is not
automatically rescored.

A valid anchor does not need to be present in every cell type. Absence from one
cell type does not cause another cell type's catalogue or cache to be used.

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
  original core set for that cell type;
- `any_direct`: match any directly associated gene in the original core set for
  that cell type.

## Mapping scope

The function includes only direct equivalents sharing a KEGG reaction ID,
Reactome reaction ID, or master Rhea ID with an anchor.

It does not perform subsystem, transitive, metabolite-neighbour, or one-hop
expansion. For cell type `t`, a mapped target is scoreable only when it is
present in every required medium-specific union GEM for `t`.

The availability rule is therefore:

```text
intersection across media within one cell type
```

It is not:

```text
intersection or union across different cell types
```

Reactions already scored as original Stage 5 cores for the same cell type are
not recomputed.

## Regulatory evidence used by the second pass

The second pass uses:

```r
step4$reaction_expression
```

For a condition-aware run this is the primary fixed-dictionary condition route.

```r
step4$reaction_expression
```

The aligned Layer 1 matrix may retain a common reaction-row universe in memory,
but each LP extracts only the reactions present in the matching cell-type union
GEM. A model for one cell type is never assigned to a metacell from another cell
type.

The targeted reactions therefore use the same condition-specific regulatory
evidence as the original Stage 5 ranking. Targets without significant estimable
condition-GRN edges retain the same RNA-only neutral fallback used in Stage 4.

## Structural reuse

For each targeted cell type and medium, the following are fixed and reused:

- the exact cell-type/medium union-GEM file and checksum;
- that cell type's biological reaction union;
- reaction order and medium-specific bounds;
- FASTCORE support previously selected inside that cell-type union GEM;
- target direction and `omega`;
- cell-type- and medium-specific directional `vmax`;
- metacell order and condition/cell-type metadata;
- Layer 1 GPR aggregation settings.

The cache contract requires:

```text
build_strategy = celltype_medium_union_gem
completion_stage = celltype_specific_fastcore_after_condition_module_union
shared_across_conditions = TRUE
shared_across_cell_types = FALSE
structural_scope = cell_type_x_medium
```

The targeted pass does not alter the original Stage 5 cache, does not merge
cached models, and does not rerun FASTCORE.

A targeted request may legitimately involve only a subset of the cell types in
Stage 4. Every reused cache cell type must exist in Stage 4, but unrelated cell
types need not be included in the targeted result.

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

Independent tasks are formed from a reused cell-type model and matching
metacells. Each worker forces numerical libraries and HiGHS to one internal
thread. Reduce the number of workers when the cached GEM or target set makes
per-worker memory the limiting resource.

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

Relation-level provenance includes cell type:

```r
targeted$expanded_reaction_catalog[, c(
  "cell_type",
  "anchor_core_reaction_id",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "present_in_merged_catalogue",
  "available_in_all_cached_union_gems"
)]
```

Target-level output also retains cell type:

```r
targeted$expanded_scoring_targets[, c(
  "cell_type",
  "reaction_id",
  "anchor_core_reaction_ids",
  "expansion_types",
  "source_annotations"
)]
```

Check exact structural reuse with:

```r
targeted$microcompass$model_cache_summary[, c(
  "cell_type",
  "medium_scenario",
  "file",
  "file_checksum",
  "reused_without_rebuilding"
)]

targeted$microcompass$params[c(
  "structural_model_reused_exactly",
  "fastcore_rerun",
  "model_rebuild",
  "structural_scope",
  "shared_across_cell_types",
  "scoring_time_limit"
)]
```

The persistent catalogue is written to
`merged_meta_module_catalogue_membership.tsv.gz`. Model cache checksums are
retained to verify exact cell-type structural reuse.

## When to rerun earlier stages

Changing only the selected anchors does not require rerunning Stages 1–5.

Changing Stage 4 reaction evidence invalidates the targeted pass. Changing a
Stage 5 medium, GEM, cell-type biological catalogue, FASTCORE completion,
bounds, direction, or target-flux control requires regenerating the affected
cell-type/medium cache before targeted scoring.

Changing an unrelated cell type does not require rebuilding cached models for a
cell type whose Stage 3 catalogue and Stage 5 structural inputs are unchanged.
