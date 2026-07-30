# RegCompassR public API

This page lists exported functions and their main purpose. See
[Workflow](workflow.md) for stage responsibilities and
[Mathematical model](mathematical-model.md) for equations.

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass_one_shot()` | Prepare defaults and run all stages. |
| `rc_run_regcompass()` | Run all stages with explicit GEM, medium, and worker settings. |

Common argument groups:

- `pando_args`: GRN fitting and filtering;
- `metacell_args`: reductions, dimensions, fixed gamma, seed, and thresholds;
- `meta_module_args`: reaction catalogue expansion inputs;
- `layer1_args`: comparison support and GPR aggregation;
- `layer2_args`: direction, solver, and model completion;
- `upstream_workers`, `layer2_workers`: stage worker budgets.

## Restartable stages

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | Pando `ConditionGRNFit v5` contracts and absolute condition coefficient tables |
| 2 | `rc_regcompass_step_metacells()` | fixed-γ condition × broad-cell-type metacells and cache provenance |
| 3 | `rc_regcompass_step_meta_modules()` | supported genes, core reactions, and merged reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | RNA support, structural-zero regulation, RNA-only fallback and reaction expression |
| 5 | `rc_regcompass_step_layer2()` | shared medium-specific models and directional scores |
| 6 | `rc_regcompass_step_results()` | annotations, rankings, evidence, and metabolic comparisons |

## Stage 1 arguments

```r
pando_args = list(
  min_cells = 100L,
  min_abs_estimate = 0,
  min_model_rsq = 0.1,
  pando_infer_args = list(
    candidate_screen = "motif_domain",
    condition_mix = 0.5,
    condition_weight = "equal",
    outer_nfolds = 5L,
    inner_nfolds = 5L,
    scale = TRUE
  )
)
```

Nested argument routing:

| Bundle | Target function |
|---|---|
| `pando_initiate_args` | `Pando::initiate_grn()` |
| `pando_motif_args` | `Pando::find_motifs()` |
| `pando_infer_args` | `Pando::infer_condition_grn()` |

RegCompass controls the object, assays, genome, metadata columns, target genes,
network name, minimum condition size, error policy, and `BPPARAM`.

Condition coefficients are absolute effects on one shared equal-condition
coordinate. No baseline-condition coefficient contrast is returned or consumed.

Important outputs:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$grn_result$normalization_policy
```

## Stage 2 arguments

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  seed = 12345L,
  min_cells_per_stratum = 100L,
  min_metacell_size = 20L,
  min_metacells_per_stratum = 2L
)
```

The same `gamma = 30L` is passed to every condition × broad-cell-type stratum.
Depth does not modify gamma. Top-1% RNA/ATAC depth cells are recorded only for
QC and never reject a metacell.

Cache validation includes cells, assays, reductions, dimensions, seed, gamma,
and thresholds.

## Stage 4 arguments

```r
rc_regcompass_step_layer1(
  grn,
  metacells,
  meta_modules,
  gem,
  outdir,
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 1,
  gpr_and_method = "min"
)
```

`projection_component` must be `"condition"`. `comparison_support` accepts
`"auto"`, `"pairwise_common"`, or `"global_common"`.
`regulatory_alpha = 1` is fixed; other values are rejected.
`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`.

The RNA Gamma–Poisson prior is estimated separately by broad cell type. A
non-estimable Pando edge contributes a structural zero in the main analysis. A
non-finite target modifier uses a neutral value and therefore equals RNA-only
support for that gene–metacell entry.

COMPASS GPR semantics are used: OR isozyme branches are summed while unavailable
branches are ignored. Missing final reaction expression is set to zero before
penalty conversion and receives the maximum expression-linked penalty.

## GEM and medium

| Function | Purpose |
|---|---|
| `rc_prepare_gem()` | Load and validate a pinned Human-GEM or Mouse-GEM. |
| `rc_make_medium_scenarios()` | Build preset or custom exchange bounds. |
| `rc_build_reaction_annotations()` | Build reaction names, formulas, GPRs, and cross-references. |
| `rc_attach_reaction_annotations()` | Attach annotations to an existing result. |

See [Medium presets](medium-presets.md).

## Condition analysis

| Function | Purpose |
|---|---|
| `rc_test_condition_reactions()` | Compare fixed reaction-direction targets across conditions. |
| `rc_report_condition_directions()` | Summarize forward and reverse targets without adding them. |
| `rc_plot_condition_reaction()` | Plot one reaction direction across conditions. |
| `rc_plot_condition_gene_reactions()` | Select and plot reaction targets associated with genes. |
| `rc_select_gene_reactions()` | Select scored reactions by metabolic gene. |

Metacell tests describe within-dataset separation and are not donor-level
inference.

## Targeted scoring

`rc_regcompass_step_target_union()` scores direct KEGG, Reactome, or master-Rhea
equivalents of selected anchors using the cached Stage 5 model. It does not
rerun FASTCORE.

## Parallel execution

- complete workflows: use `upstream_workers` and `layer2_workers`;
- stepwise functions: pass a `BiocParallelParam` through `BPPARAM`;
- `BPPARAM = TRUE` is invalid;
- do not set nested `pando_infer_args$parallel` in complete-workflow examples.
