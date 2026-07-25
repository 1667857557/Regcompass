# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.4.

## Workflow

```text
condition × cell type cells
→ Pando GRNs and multimodal metacells
→ complete-GPR reaction meta-modules
→ integrated RNA+ATAC reaction support
→ medium-specific model and global FASTCORE completion
→ directional LP scoring and condition contrasts
```

## Prepare the object and model

```r
library(RegCompassR)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motif2tf, package = "Pando")

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = "high_glucose",
  species = "human"
)
```

The Seurat object must contain normalized RNA and ATAC assays and the metadata columns supplied below. Pando is fitted separately for each `condition × cell type` group.

Available medium presets include physiological plasma, RPMI-1640, high-glucose DMEM, glucose/lactate/glutamine sensitivity scenarios, technical exchange baselines, and custom media. See [medium presets](medium-presets.md) for the complete list and assumptions.

## Run the complete workflow

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motif2tf,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  gem = gem,
  species = "human",
  medium_scenarios = medium_scenarios,
  sample_col = NULL,
  condition_col = "Group",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  metacell_args = list(
    gamma = 30,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  pando_args = list(
    min_cells = 100,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    max_targets_per_tf = 200,
    expansion_mode = "ordered_once",
    regulatory_alpha = 1,
    tau = 0.20
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 600,
    model_params = list(
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  ),
  upstream_workers = 6,
  layer2_workers = 30,
  progress = TRUE
)
```

`layer2_args$model_params` controls the medium-specific FASTCORE completion. One cached structural model is shared by all conditions analysed under the same medium.

## Inspect the main outputs

```r
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```

Within one medium scenario, condition differences arise from the RNA+ATAC penalty matrix rather than condition-specific network reconstruction.
