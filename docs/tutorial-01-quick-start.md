# Tutorial 1: minimal one-shot run

Use this tutorial for a complete analysis from a paired-cell RNA+ATAC Seurat
object. Mathematical definitions are in
[Mathematical model](mathematical-model.md).

## Required object state

The object must contain:

- paired RNA and ATAC assays for the same cells;
- condition and broad-cell-type metadata;
- RNA PCA and ATAC LSI reductions;
- genome-compatible ATAC peak coordinates.

```r
stopifnot(
  all(c("Group", "cell_type") %in% colnames(A@meta.data)),
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions)
)
```

RegCompass selects GEM GPR genes present in the RNA assay. Do not override
`genes` inside `pando_infer_args`.

## Prepare GEM and medium

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
```

Use `physiologic` for the species-specific baseline. Culture-medium and nutrient
challenge presets are intended for sensitivity analysis. Use `custom` for the
actual experimental environment. See [Medium presets](medium-presets.md).

## Run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100L,
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
    min_cells_per_stratum = 500L,
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

Do not add `parallel` or `BPPARAM` inside `pando_infer_args`; the one-shot runner
controls Stage 1 workers.

## Key parameters

- `candidate_screen = "motif_domain"`: canonical candidate policy.
- `reference_condition`: controls interpretation contrasts, not the primary penalty.
- `gamma`: approximate cells per metacell.
- `regulatory_alpha`: strength of regulatory modification of RNA support.
- `gpr_and_method = "min"`: limiting-subunit GPR rule.
- `target_direction = "both"`: score forward and reverse directions separately.

## Inspect outputs

```r
result$grn$condition_fit_status
result$grn$condition_grn_fits
result$grn$tf_peak_gene_condition
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$microcompass$model_cache_summary
result$reaction_ranking
result$condition_contrast
```

The primary reaction comparison uses outer-heldout common-support regulatory
scores. Reference-condition coefficient effects and condition-full projections
are interpretation or sensitivity outputs.

## Mouse input

Mouse analyses require build-matched regulatory regions:

```r
library(BSgenome.Mmusculus.UCSC.mm10)

mouse_regions <- readRDS("mm10_regulatory_regions.rds")

mouse_result <- rc_run_regcompass_one_shot(
  object = A_mouse,
  outdir = "RegCompass_mouse",
  genome = BSgenome.Mmusculus.UCSC.mm10,
  species = "mouse",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    pando_initiate_args = list(regions = mouse_regions),
    pando_infer_args = list(reference_condition = "Control")
  )
)
```

The region build must match the ATAC coordinates and motif-scanning genome.

Use [Tutorial 2](tutorial-02-stepwise-audit.md) for restartable stages. Public
API: [functions.md](functions.md).
