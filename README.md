# RegCompassR

RegCompassR connects paired single-cell RNA+ATAC regulatory evidence to
metacell-level GEM and COMPASS-like reaction scoring.

## Canonical workflow

```r
library(RegCompassR)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L
  ),
  layer1_args = list(
    projection_component = "condition",
    comparison_support = "auto",
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),
  medium_scenarios = medium_scenarios,
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  )
)
```

With at least two conditions, Stage 1 uses `condition_grn`; otherwise it uses
`standard_pando` through `Pando::infer_grn()` and calculates No condition
coefficients.

A canonical run may explicitly omit condition metadata:

```r
single_result <- rc_run_regcompass(
  object = A,
  gem = gem,
  outdir = "RegCompass_single",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = NULL,
  celltype_col = "cell_type"
)
```

In this route, `single_result$reaction_ranking` remains available and
`single_result$condition_contrast` is empty.

For condition-aware analysis:

- each broad cell type has one shared TF–peak–target candidate supergraph;
- the primary regulatory signal is `condition_full_oof`;
- jointly estimable edges form the common-support component;
- an edge non-estimable in one or both conditions contributes a projectable
  structural zero in each affected condition;
- exact-zero predictors remain represented without receiving fitted
  coefficients.

Stage 2 calls `SCimplify_by_graph_group_from_embedding()` with
`cell.graph.group = cell_type` and `cell.split.condition = condition`. The graph
scope is `one_independent_graph_per_cell_type`; all conditions are joint within
that graph (`all_conditions_joint_within_cell_type_graph`) and
`temporary_combined_stratum = FALSE`.

The primary metabolic ranking uses the condition-full penalty. Common-support
and RNA-only penalties are retained as decomposition/control outputs. The
canonical schema does not calculate depth matching, common-depth restriction,
alpha sensitivity, zero-support sensitivity, or link-saturation propagation.

## Optional targeted reaction remapping

After a completed `model_mode = "meta_module_gem"` run,
`rc_regcompass_step_target_union()` can score direct KEGG-, Reactome-, or
master-Rhea-linked non-core reactions. It reuses the exact cached Stage 5 union
GEM and global FASTCORE completion. It is a retained optional second LP pass,
not a comparability guardrail, and uses the canonical Layer 1
`reaction_expression` input, which is
`reaction_expression_condition_full_oof` in condition mode.

## Documentation

- [Tutorial 1: one-shot workflow](docs/tutorial-01-quick-start.md)
- [Tutorial 2: stepwise workflow](docs/tutorial-02-stepwise-audit.md)
- [Tutorial 3: mathematical model](docs/tutorial-03-mathematical-model.md)
- [Tutorial 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Tutorial 5: condition comparison](docs/tutorial-05-condition-differential-analysis.md)
- [Public functions](docs/functions.md)
- [Stage schemas](docs/stage-interface-contracts.md)

Metacells are valid within-dataset statistical units. Their P values quantify
condition-associated metacell separation and are not donor-level biological
replicate inference.
