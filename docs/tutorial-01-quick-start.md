# Tutorial 1: quick start

## Input

Use a paired-cell Seurat object with RNA and ATAC count assays, broad cell-type metadata, RNA PCA, ATAC LSI and genome-compatible peak coordinates. Stage 1 uses `pando_args$min_cells = 500L` by default for each retained condition-by-cell-type stratum, but this threshold is user-configurable. Stage 2 likewise defaults `metacell_args$min_cells_per_stratum` to `500L` and allows an explicit override.

## GEM and medium

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
```

Built-in biological scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine` and `custom`. See [medium-presets.md](medium-presets.md) for custom reaction bounds, metabolite compositions and publication provenance.

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
    min_cells = 500L,
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
    min_cells_per_stratum = 500L,
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
  workers = 10L
)
```

The two `500L` cell-count settings are defaults, not fixed constraints. If a dataset requires a lower or higher threshold, set `pando_args$min_cells` and/or `metacell_args$min_cells_per_stratum` explicitly; RegCompass preserves those supplied values.

`workers` is the only workflow-level parallel setting. Its default is `10L` and it may be changed, for example `workers = 60L`. RegCompass automatically selects `SnowParam(type = "SOCK")` on Windows and `MulticoreParam` on Linux/macOS. The effective cap is `min(workers, max(1, detected logical CPUs - 2))`, and each individual Pando/CORDA2/LP dispatch shrinks further to its own independent task count.

For cell types with at least two retained conditions, RegCompass first parallelizes pooled-background and condition × cell-type candidate discovery, freezes one exact edge dictionary per cell type after a strict barrier, and then parallelizes condition × cell-type fixed-dictionary GLMs. This tutorial uses `tf_cor = 0.1` and `peak_cor = 0.05` for that candidate screen. Each atomic Pando task disables nested target-level parallelism. If only one effective condition is retained, standard Pando is parallelized across broad cell types and uses the same fixed `tf_cor` and `peak_cor` values supplied through `pando_infer_args`; RegCompass does not increase either threshold as a function of cell count.

`meta_module_gem` uses original MATLAB CORDA2 by default. Set `model_params$model_completion = "fastcore"` only for the supplementary FASTCORE route. Use `model_mode = "full_gem"` for supplementary complete-network COMPASS-style scoring.

CORDA2 reconstruction intentionally runs without a structural time limit. Do not supply `model_params$completion_time_limit` for the default CORDA2 route; that control is reserved for supplementary non-CORDA2 completion such as FASTCORE.

CORDA2 receives the complete medium-constrained parent without FASTCC pre-pruning. Retained reactions recover their parent directional bounds, including positive lower bounds. Layer 2 uses the COMPASS cost scale; missing expression and structural roles receive cost `1`. Step 1, Step 2.1, Step 2.2 and Step 3 remain strict mathematical barriers; directional targets within each step are parallelized up to the protected worker cap, the step pool is then released, and the next step starts with a fresh pool.

## Main outputs

```r
result$grn$cell_type_analysis_mode
result$layer1$gene_regulatory_modifier
result$microcompass$penalty
result$reaction_ranking
result$condition_contrast
```

Detailed stage parameters are in [tutorial-02-stepwise-audit.md](tutorial-02-stepwise-audit.md). Mathematical definitions are in [mathematical-model.md](mathematical-model.md).
