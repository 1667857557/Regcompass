# Tutorial 1: quick start

## Input requirements

Use a paired-cell Seurat object with RNA and ATAC counts, broad cell-type metadata, RNA PCA, ATAC LSI, and genome-compatible peak coordinates.

Stage 1 uses a fixed minimum of 300 paired cells. For each cell type:

- at least two retained conditions: common-dictionary condition GRN;
- one retained condition: standard Pando;
- no retained condition stratum: excluded.

The same `pando_infer_args` list can be supplied in all cases. Stage 1 automatically disables condition-only controls before standard Pando and disables standard-model controls for the fixed common-dictionary condition model. With multiple retained cell types, `upstream_workers` distributes independent cell-type GRN jobs and prevents nested worker pools.

RegCompass defaults both Pando candidate-discovery correlation thresholds to `0.05`. Because the current Pando implementation applies strict `>` comparisons, the effective candidate rules are:

```r
abs(tf_target_cor) > 0.05
abs(peak_target_cor) > 0.05
```

Candidate discovery and final penalty entry are separate operations. A fitted TF-peak-target edge contributes to the RegCompass regulatory penalty only when all four conditions hold:

```r
estimable == TRUE &
  padj < 0.05 &
  abs(corr) >= 0.05 &
  abs(estimate) >= 0.05
```

The two final absolute-value gates are inclusive. Positive and negative correlations and coefficients are retained when their absolute magnitudes pass the threshold. Edges that fail any final gate remain in the complete coefficient table but receive `penalty_effect = 0` and do not contribute to downstream regulatory projection. For standard Pando, `corr` is the coefficient-table TF-target correlation. For a common-dictionary condition fit without coefficient-level `corr`, RegCompass uses the frozen dictionary's audited `max_abs_tf_target_cor` and records the origin in `corr_source`.

## Install current companion repositories

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Prepare GEM and medium

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

Supported built-in scenarios include:

| Scenario | Use |
|---|---|
| `normal_human_plasma` | Human plasma-like background |
| `mouse_plasma` | Mouse plasma/interstitial-fluid background |
| `high_glucose` | Plasma-like background with increased glucose |
| `low_glucose` | Plasma-like background with reduced glucose |
| `high_lactate` | Plasma-like background with increased lactate |
| `low_lactate` | Plasma-like background with reduced lactate |
| `low_glutamine` | Plasma-like background with reduced glutamine |
| `custom` | User-supplied exchange bounds or metabolite availability |

Custom reaction bounds:

```r
custom_medium <- data.frame(
  medium_scenario_id = "measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

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
      tf_cor = 0.05,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
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

## Main outputs

```r
result$grn$cell_type_analysis_mode
result$grn$condition_fit_status
result$grn$pando_execution_plan
result$grn$pando_infer_argument_routing
result$layer1$gene_regulatory_modifier
result$microcompass$penalty
result$reaction_ranking
result$condition_contrast
```

The mathematical definitions are maintained only in [mathematical-model.md](mathematical-model.md).
