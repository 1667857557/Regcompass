# RegCompassR

RegCompassR 2.1.0 integrates paired single-cell RNA+ATAC regulatory evidence
with a shared metabolic model for condition-level reaction comparison.

## Workflow

```text
condition-aware Pando GRNs
→ condition × broad-cell-type metacells
→ supported metabolic genes and reactions
→ regulatory modification of RNA support
→ shared medium-specific metabolic model
→ directional reaction scores and condition tests
```

Pando is the GRN estimator. RegCompass uses Pando `ConditionGRNFit v5`
outer-heldout common-support projections for the primary penalty and retains
reference-condition effects for interpretation only.

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
      reference_condition = "Control",
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
    regulatory_alpha = 0.5,
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

## Important choices

### Candidate screening

RegCompass requires `candidate_screen = "motif_domain"`. Alternative Pando
candidate-screen modes are not accepted by the current Stage 1 wrapper.

### GPR AND aggregation

- `min`: default limiting-subunit rule;
- `median`, `mean`: sensitivity options.

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
result$microcompass$model_cache_summary
result$reaction_ranking
result$condition_contrast
result$reaction_comparison_by_metacell
```

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
- `rc_regcompass_step_metacells()` — condition × broad-cell-type metacells.
- `rc_regcompass_step_meta_modules()` — supported genes and reaction catalogue.
- `rc_regcompass_step_layer1()` — gene support and reaction penalties.
- `rc_regcompass_step_layer2()` — shared models and directional scores.
- `rc_regcompass_step_results()` — rankings, annotations, and contrasts.

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
