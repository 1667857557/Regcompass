# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.4.

## Current architecture

```text
condition × cell type single cells
→ Pando GRN
→ complete-GPR core reactions
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ biological meta-modules
→ deduplicated merged meta-module catalogue
→ Layer 1 RNA+ATAC reaction support
→ one medium-specific union GEM
→ one global FASTCORE completion
→ directional COMPASS-like LP scoring
```

Stage 3 does **not** run FASTCORE and does **not** create a GEM. The phrase **union GEM** is reserved for the medium-constrained model created in Stage 5 after all biological meta-modules have been merged.

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

Do not supply `local_fastcore` or `local_fastcore_args`. Those interfaces were removed. FASTCORE is configured only through `layer2_args$model_params` and is applied once to each medium-specific union GEM.

## Inspect the main outputs

```r
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```

`merged_grn_meta_modules` is a reaction catalogue, not a GEM. To inspect an actual union GEM:

```r
union_summary <- result$microcompass$model_cache_summary
union_gem <- readRDS(union_summary$file[[1]])

union_gem$is_union_gem
union_gem$union_gem_medium_scenario
union_gem$build_params[c(
  "n_biological_reactions",
  "n_fastcore_support_reactions",
  "completion_stage"
)]
```

All conditions and metacells evaluated under the same medium scenario use the same cached union GEM. Their differences arise from the RNA+ATAC penalty matrix, not condition-specific network structures.
