# Tutorial 1: quick start

## Input

Use a paired-cell Seurat object with RNA and ATAC count assays, broad cell-type metadata, RNA PCA, ATAC LSI and genome-compatible peak coordinates. Stage 1 requires at least 300 paired cells per retained condition-by-cell-type stratum.

## GEM, medium and worker budget

```r
library(RegCompassR)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)

workers <- 10L
```

Built-in biological scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine` and `custom`. See [medium-presets.md](medium-presets.md) for custom reaction bounds, metabolite compositions and publication provenance.

`workers` is the only workflow-wide parallel setting. Ten is the default request and can be adjusted. RegCompass automatically chooses SOCK workers on Windows and multicore workers on Linux/macOS, reserves two detected CPUs, and caps the resolved budget at `max(1, available CPUs - 2)`. Each operation then uses no more than its independent task count.

## Run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    min_cells_per_stratum = 300L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  ),
  layer1_args = list(
    gpr_and_method = "min",
    gene_half_saturation = 1
  ),
  medium_scenarios = medium_scenarios,
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      strict = TRUE,
      corda2_args = list(
        MCxNCthresh = 2,
        constraint = 1,
        constrainby = "val",
        om = 1e4,
        ci = 0.01
      ),
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20
    )
  ),
  workers = workers
)
```

For cell types with at least two retained conditions, the tutorial uses `tf_cor = 0.1` for the common-dictionary condition-GRN route. If only one effective condition is retained, standard Pando treats the requested `tf_cor` as an effect-size floor and automatically raises it when necessary to the two-sided Pearson-correlation critical value for that cell type's actual cell count (`alpha = 0.05`). This sample-size-aware gate prevents small cell groups from receiving a more permissive TF-correlation screen solely because their sampling variance is larger.

Stage 1 parallelizes independent broad cell types and therefore uses only `min(number of cell-type jobs, resolved workers)` workers. CORDA2 and large directional LP target sets can use the complete resolved budget when enough targets are present.

`meta_module_gem` uses original MATLAB CORDA2 by default. Set `model_params$model_completion = "fastcore"` only for the supplementary FASTCORE route. Use `model_mode = "full_gem"` for supplementary complete-network COMPASS-style scoring.

CORDA2 reconstruction intentionally runs without a structural time limit. Do not supply `model_params$completion_time_limit` for the default CORDA2 route; that control is reserved for supplementary non-CORDA2 completion such as FASTCORE.

CORDA2 receives the complete medium-constrained parent without FASTCC pre-pruning. Retained reactions recover their parent directional bounds, including positive lower bounds. Layer 2 uses the COMPASS cost scale; missing expression and structural roles receive cost `1`.

Within each CORDA2 mathematical step, directional targets are parallelized behind strict barriers. The console and Layer 2 progress files report completed and remaining target counts. The step-local worker pool is released and full garbage collection runs before the next CORDA2 step starts.

## Main outputs

```r
result$grn$cell_type_analysis_mode
result$layer1$gene_regulatory_modifier
result$microcompass$penalty
result$reaction_ranking
result$condition_contrast
```

Detailed stage parameters are in [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). Mathematical definitions are in [mathematical-model.md](mathematical-model.md).
