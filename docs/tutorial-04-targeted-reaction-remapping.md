# Tutorial Level 4: targeted reaction remapping on the cached union GEM

**Previous:** [Tutorial 2](tutorial-02-stepwise-audit.md) creates the required `step3`, `step4`, and `step5` objects. [Tutorial 3](tutorial-03-advanced-restart.md) explains when those objects can be reused.

**This tutorial:** scores selected non-core reactions that are directly equivalent to user-specified reaction anchors. It reuses the exact Stage 5 union GEM and does not rerun FASTCORE.

**Next:** [Tutorial 5](tutorial-05-condition-differential-analysis.md) applies the same condition-analysis logic to either original Stage 5 targets or this second-pass target set.

## 1. Load matching stage checkpoints

```r
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

They must come from the same validated workflow and GEM fingerprint. `step5` must use:

```r
stopifnot(identical(step5$model_mode, "meta_module_gem"))
step5$model_cache_summary
```

The cached model already contains the final medium constraints, global FASTCORE support, reaction order, `S`, `lb`, and `ub`.

## 2. Select anchors by reaction ID

```r
targeted <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_targeted/reaction_anchors",
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

The compatibility argument `core_reaction_ids` accepts either an original core or another reaction in the supplied GEM. A reaction-ID anchor is used to discover direct database equivalents; it is not automatically added as a new score target.

## 3. Select original core anchors by gene

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

`gene_match = "complete_gpr"` requires at least one complete original GPR branch. `"any_direct"` is a less stringent direct-gene match within the original core set and should be treated as an exploratory anchor rule.

## 4. Exact mapping scope

The second pass includes only non-core reactions sharing at least one direct identifier with an anchor:

```text
KEGG reaction ID
OR Reactome reaction ID
OR master Rhea ID
```

It does not perform:

```text
subsystem expansion
transitive cross-reference expansion
metabolite-neighbour expansion
one-hop network expansion
new FASTCORE completion
```

A mapped reaction is scored only when it already exists in every cached union GEM required by the requested media. Therefore the comparison keeps the same structural model used in the original analysis.

## 5. Inspect relation-level provenance

```r
targeted$selected_anchor_reactions
targeted$selected_core_reactions
targeted$selected_noncore_reactions
targeted$expanded_reaction_catalog
targeted$expanded_scoring_targets
```

```r
relation_provenance <- targeted$expanded_reaction_catalog[, intersect(c(
  "anchor_reaction_id",
  "anchor_is_original_core",
  "reaction_id",
  "expansion_type",
  "source_annotation",
  "present_in_merged_catalogue",
  "merged_catalogue_inclusion_stage",
  "available_in_all_cached_union_gems"
), colnames(targeted$expanded_reaction_catalog)), drop = FALSE]

relation_provenance
```

Reaction-level scoring targets:

```r
targeted$expanded_scoring_targets[, intersect(c(
  "reaction_id",
  "anchor_reaction_ids",
  "expansion_types",
  "source_annotations",
  "merged_catalogue_inclusion_stage"
), colnames(targeted$expanded_scoring_targets)), drop = FALSE]
```

## 6. Verify that the structural model was reused

```r
targeted$microcompass$model_cache_summary
targeted$microcompass$params[c(
  "structural_model_reused_exactly",
  "fastcore_rerun",
  "model_rebuild",
  "scoring_time_limit"
)]
```

Expected contract:

```text
structural_model_reused_exactly = TRUE
fastcore_rerun = FALSE
model_rebuild = FALSE
scoring_time_limit = absent/NULL
```

The original `completion_time_limit` applied only during Stage 5 model construction. It is neither reused nor configurable in this scoring-only second pass.

## 7. Compare original and second-pass result families

Keep source labels distinct:

```r
original_report <- rc_report_condition_directions(
  step5,
  source_label = "original_stage5_targets"
)

second_pass_report <- rc_report_condition_directions(
  targeted,
  source_label = "target_union_direct_equivalents"
)
```

Do not silently combine adjusted P values from original and second-pass target families. Either report them separately or explicitly define and recompute a combined multiple-testing family.

## 8. Persistent outputs

The target-union stage writes its own restart object and mapping tables under the supplied `outdir`. The persistent merged catalogue provenance includes:

```text
merged_meta_module_catalogue_membership.tsv.gz
```

The original compact `result` remains unchanged. The second pass is a separate scoring object, not a mutation of Stage 3 core definitions.

## 9. Handoff to Tutorial 5

In [Tutorial 5](tutorial-05-condition-differential-analysis.md), use:

- `result` for original Stage 5 core/union targets;
- `targeted` for direct-equivalent second-pass targets;
- a distinct `source_label` for each;
- the same medium, cell type, condition pair, and direction when comparing reports.
