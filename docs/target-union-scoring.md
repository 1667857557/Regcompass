# Direct database-linked non-core scoring in cell-type union GEMs

`rc_regcompass_step_target_union()` is an optional second-pass analysis after a
completed stepwise Layer 2 run with `model_mode = "meta_module_gem"`. It reuses
the exact cell-type-specific structural models built during the primary run and
does not create a new cross-cell-type reaction universe.

## Structural source

The function reuses the final union-GEM files recorded by:

```r
layer2$model_cache_summary
```

Each cache row identifies one `cell_type × medium_scenario` model and records:

- cell type;
- medium scenario;
- model file and checksum;
- build strategy;
- FASTCORE completion stage.

Each cached model contains only:

- biological meta-module reactions unioned across conditions of that cell type;
- medium-specific bounds for that cell type;
- FASTCORE support selected independently for that cell-type model;
- the original direction-specific core targets for that cell type.

The Stage 3 object supplies cell-type-scoped core and membership tables through:

```r
meta_modules$merged_modules$merged_core_reactions
meta_modules$merged_modules$merged_reaction_membership
meta_modules$merged_modules$cell_type_catalogues
```

These Stage 3 tables are biological catalogues, not GEMs. Every core and
membership row carries the workflow cell-type column.

## Evidence route

The second pass obtains reaction expression from:

```r
layer1$reaction_expression
```

For condition-aware runs, this is the canonical alias of:

```r
layer1$reaction_expression_condition_full_oof
```

The aligned expression matrix may contain a common reaction-row universe in
memory, but each LP uses only the reactions present in the matching cell-type
GEM. A model for one cell type is never assigned to a metacell from another cell
type.

## Mapping rule

Reaction IDs supplied through `core_reaction_ids` are mapping anchors. An anchor
is evaluated separately in every cell type in which it belongs to the original
core or is present in the cached structural model.

Direct candidate targets are reactions sharing at least one:

- KEGG reaction identifier;
- Reactome reaction identifier;
- master Rhea identifier.

`core_genes` resolves anchors within each cell type's original complete-GPR core
set. No subsystem expansion, transitive expansion, metabolite-neighbour
expansion, model reconstruction, or FASTCORE completion is performed.

## Cell-type availability rule

A mapped non-core reaction is scoreable for cell type `c` only when it is present
in every required medium-specific union GEM for `c`. Availability is intersected
across media **within one cell type**. It is never intersected across different
cell types.

Before scoring, RegCompass verifies:

```text
cell_type matches the cache row and model provenance
file_checksum matches
build_strategy = celltype_medium_union_gem
completion_stage = celltype_specific_fastcore_after_condition_module_union
model$is_union_gem = TRUE
model$union_gem_scope =
  one_cell_type_one_medium_shared_across_conditions_and_matching_metacells
model$union_gem_medium_scenario matches the cache row
```

FASTCORE support reactions may be targeted when they are present in all reused
media for the same cell type, even when absent from that cell type's Stage 3
biological catalogue. Original Layer 2 core targets are not recomputed.

## Time-limit policy

`layer2_args$model_params$completion_time_limit` belongs to the original Stage 5
construction of each cell-type union GEM. The targeted pass reuses those files
without model reconstruction or FASTCORE reruns, so it has no construction
`time_limit` parameter.

## Example

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted",
  core_reaction_ids = c(
    "MAR04381",
    "MAR04379",
    "MAR04391"
  ),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
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

The relation and target tables include the workflow cell-type column together
with:

```text
anchor_core_reaction_id
reaction_id
expansion_type
source_annotation
present_in_merged_catalogue
merged_catalogue_is_core
merged_catalogue_inclusion_stage
available_in_all_cached_union_gems
```

The scoring result records exact cache reuse and cell-type scope through:

```r
targeted$microcompass$params[c(
  "structural_model_reused_exactly",
  "fastcore_rerun",
  "model_rebuild",
  "structural_scope",
  "shared_across_cell_types",
  "scoring_time_limit"
)]
```

Tutorial: [targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md).
