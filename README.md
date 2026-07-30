# RegCompassR

RegCompassR 2.1.0 integrates paired single-cell RNA+ATAC regulatory evidence
with a shared metabolic model for condition-level reaction comparison.

## Workflow

```text
absolute condition-aware Pando GRNs
→ condition × broad-cell-type fixed-γ metacells
→ cell-type-specific RNA latent-expression priors
→ structural-zero regulatory projection + RNA-only fallback
→ COMPASS GPR reaction expression and penalties
→ shared medium-specific metabolic model
→ directional reaction scores and condition tests
```

Pando is the GRN estimator. RegCompass uses Pando `ConditionGRNFit v5`
outer-heldout condition projections for the primary penalty. Conditions are
represented by absolute coefficients on one equal-condition coordinate; no
reference-condition coefficient or stored GRN contrast enters RegCompass.

Detailed equations: [Mathematical model](docs/mathematical-model.md).
Stage descriptions: [Workflow](docs/workflow.md).

## Installation

### Validated Seurat v4 profile

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@Supercell2")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

A coherent SeuratObject/Seurat 5.x stack with Signac 1.12–1.x is also accepted.
Signac 2.x `ChromatinAssay5` is not supported. See
[Seurat compatibility](docs/seurat-compatibility.md).

## Required input

A paired-cell Seurat object containing:

- RNA and ATAC assays for the same cells;
- condition and broad-cell-type metadata;
- RNA PCA and ATAC LSI reductions for metacell construction;
- genome-compatible peak coordinates.

Mouse analyses must supply build-matched regulatory regions through
`pando_initiate_args$regions`.

## Minimal one-shot run

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = "physiologic",
  species = "human"
)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  cell_type = "T_cell",
  pando_args = list(
    min_cells = 300L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      condition_mix = 0.5,
      condition_weight = "equal",
      outer_nfolds = 5L,
      inner_nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 300L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
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
    solver = "highs",
    model_params = list(
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  ),
  upstream_workers = 6L,
  layer2_workers = 30L
)
```

The one-shot workflow controls Stage 1 parallelism through `upstream_workers`.
Do not set `parallel` or `BPPARAM` inside `pando_infer_args`.

## Canonical analysis policies

### Metacells

Every condition × broad-cell-type stratum uses the same `gamma`, default 30.
`gamma` is not adjusted by RNA/ATAC depth. Cells above the stratum 99th
percentile of RNA or ATAC depth are counted for diagnostics only; no metacell is
rejected because it contains more than one such cell.

### RNA prior

The Gamma–Poisson empirical-Bayes prior is estimated independently for every
broad cell type. Conditions within the same cell type share that prior, while
different cell types do not.

### Regulatory projection

A Pando edge that is unavailable under the requested estimability/support policy
contributes exactly zero at the single-cell projection-contribution layer. This
structural zero enters target summation, metacell averaging, GPR aggregation,
reaction expression and the main penalty path; the structural-zero mask remains
available for audit.

When a target-level Pando modifier is otherwise unavailable or non-finite,
RegCompass uses a neutral modifier (`R = 0`), making the result exactly equal to
the RNA-only support for that gene–metacell entry. The fallback is recorded in
`result$layer1$regulatory_fallback`.

`regulatory_alpha` is fixed at 1. Other values are rejected by the public Layer 1
API.

### GPR aggregation and penalty

- AND complexes use `min` by default; `median` and `mean` remain sensitivity
  options. COMPASS missing-value semantics are retained for each AND method.
- OR isozyme branches are summed while unavailable branches are ignored, matching
  COMPASS isoform summing.
- If final reaction expression is unavailable, it is set to zero before penalty
  conversion. It therefore receives the maximum expression-linked penalty of 1.
  Structural exchange/demand/sink/support reactions retain their fixed costs.

### Candidate screening

RegCompass requires `candidate_screen = "motif_domain"`. Alternative Pando
candidate-screen modes are not accepted by the current Stage 1 wrapper.

### Medium

Use `physiologic` for the species-specific baseline, a culture-medium preset for
in-vitro sensitivity analysis, or `custom` for experiment-specific bounds.
Concentrations are not measured uptake fluxes. See
[Medium presets](docs/medium-presets.md).

## Main outputs

```r
result$grn$condition_grn_fits
result$grn$condition_fit_status
result$grn$tf_peak_gene_condition
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$layer1$projection_structural_zero
result$layer1$regulatory_fallback
result$layer1$capacity_params
result$microcompass$model_cache_summary
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

`result$condition_contrast` is a downstream metabolic comparison between
conditions. It is not a Pando reference-condition coefficient contrast.

## Condition tests

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "Group",
  celltype_col = "cell_type",
  reaction_ids = "MAR06231",
  target_directions = "forward",
  medium_scenarios = "physiologic"
)

condition_stats$pairwise
condition_stats$omnibus
```

Compare the same reaction, direction, medium, and broad cell type. Metacell P
values describe within-dataset separation and are not donor-level inference.

## Restartable stages

- `rc_regcompass_step_grn()` — condition-aware Pando models.
- `rc_regcompass_step_metacells()` — fixed-γ condition × broad-cell-type metacells.
- `rc_regcompass_step_meta_modules()` — supported genes and reaction catalogue.
- `rc_regcompass_step_layer1()` — RNA support, structural-zero regulation and
  COMPASS reaction expression.
- `rc_regcompass_step_layer2()` — shared models and directional scores.
- `rc_regcompass_step_results()` — rankings, annotations, and metabolic condition
  comparisons.

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Stepwise workflow](docs/tutorial-02-stepwise-audit.md)
- [Restart and sensitivity](docs/tutorial-03-advanced-restart.md)
- [Targeted reaction scoring](docs/tutorial-04-targeted-reaction-remapping.md)
- [Condition comparison](docs/tutorial-05-condition-differential-analysis.md)
- [Workflow](docs/workflow.md)
- [Mathematical model](docs/mathematical-model.md)
- [Pando condition contract](docs/condition-comparable-grn.md)
- [Medium presets](docs/medium-presets.md)
- [Public API](docs/functions.md)
