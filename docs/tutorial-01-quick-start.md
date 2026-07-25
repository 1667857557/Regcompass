# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.4.

## Workflow

```text
condition × cell type cells
→ Pando models of Human-GEM GPR genes
→ significantly supported metabolic target genes
→ complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ integrated RNA+ATAC reaction support
→ medium-constrained model with global FASTCORE completion
→ directional LP scoring and condition contrasts
```

## Prepare the object and model

```r
library(RegCompassR)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)

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

The Seurat object must contain normalized RNA and ATAC assays and the metadata columns supplied below. Pando is fitted separately for each `condition × cell type` group. Its target list is the intersection of Human-GEM GPR genes and RNA-assay row names.

When `pfm` is omitted, RegCompass loads `data("motifs", package = "Pando")` and passes the resulting `motifs` object to `Pando::find_motifs()`. A user-supplied `pfm` overrides this default.

By default, human analyses also load the Pando data objects `phastConsElements20Mammals.UCSC.hg38` and `SCREEN.ccRE.UCSC.hg38`, take their union, and pass that `GRanges` object to `Pando::initiate_grn(regions = ...)`. Override this only through `pando_args$pando_initiate_args$regions`. Non-human analyses must provide species-appropriate regions explicitly.

Available medium presets include physiological plasma, RPMI-1640, high-glucose DMEM, glucose/lactate/glutamine sensitivity scenarios, technical exchange baselines, and custom media. See [medium presets](medium-presets.md) for the complete list and assumptions.

## Run the complete workflow

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
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
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),
  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),
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
  upstream_workers = 6,
  layer2_workers = 30,
  progress = TRUE
)
```

The Stage 1 evidence filter defines the Stage 3 gene set. A target gene is supported when at least one TF–peak–gene row passes the adjusted-P-value, absolute-estimate, and model-R² filters. Both positive and negative coefficients count as regulatory evidence.

Stage 3 always performs exactly one ordered expansion pass:

```text
core subsystem
→ KEGG/Reactome reaction equivalence
→ master-Rhea reaction equivalence
```

The retired `expansion_mode`, `max_iterations`, fixed-point, and one-hop reaction APIs have been removed.

`layer1_args$gpr_and_method` controls genes joined by a GPR AND relationship. Allowed values are `"min"`, `"median"`, and `"mean"`; the default is `"min"`. The retired Boltzmann soft-min and `tau` parameter have been removed. Isozyme OR branches are summed in the canonical Layer 1 calculation.

`layer2_args$model_params$completion_time_limit` applies only to FASTCORE union-GEM construction; scoring LPs have no time-limit parameter.

## Inspect the main outputs

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$microcompass$model_cache_summary
```
