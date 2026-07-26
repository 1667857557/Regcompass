# Tutorial Level 1: minimal one-shot run

Use this tutorial for a paired-cell RNA+ATAC Seurat object and RegCompassR 1.8.4.

## Workflow

```text
condition × cell type cells
→ Pando models of GEM GPR genes
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

The Seurat object must contain normalized RNA and ATAC assays, the requested reductions, and the metadata columns supplied below. Pando is fitted separately for each `condition × cell type` group. Its target list is the intersection of GEM GPR genes and RNA-assay row names.

When `pfm` is omitted, RegCompass loads `data("motifs", package = "Pando")` and passes the resulting `motifs` object to `Pando::find_motifs()`. A user-supplied `pfm` overrides this default.

The default Pando regions are species-specific:

```text
human: phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse: phastConsElements20Mammals.UCSC.hg38 only
```

Both objects are loaded from Pando. An explicit `pando_args$pando_initiate_args$regions` overrides the default.

Available medium presets include physiological plasma, RPMI-1640, high-glucose DMEM, glucose/lactate/glutamine sensitivity scenarios, technical exchange baselines, and custom media. See [medium presets](medium-presets.md) for the complete list and assumptions.

## Run the complete workflow

The public runner arguments and the example below follow the processing sequence: model/shared inputs → Stage 1 Pando → Stage 2 metacells → Stage 3 meta-modules → Stage 4 Layer 1 → Stage 5 Layer 2 → execution controls.

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,

  # Stage 1: Pando GRN evidence
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE,
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    )
  ),

  # Stage 2: metacells
  sample_col = NULL,
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  # Stage 4: integrated evidence and GPR aggregation
  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),

  # Stage 5: medium-specific union GEM and scoring
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

  upstream_workers = 6,
  layer2_workers = 30,
  progress = TRUE
)
```

## Pando evidence-filter parameters

A TF–peak–target row is retained only when all configured requirements are met:

| Parameter | Default | Meaning |
|---|---:|---|
| `min_cells` | `20L` | Minimum cells required for each `condition × cell type` Pando fit. The example uses `100`. |
| `padj_threshold` | `0.05` | Maximum adjusted P value for a TF–peak–target coefficient. |
| `min_abs_estimate` | `0` | Minimum absolute Pando coefficient. `0` retains every finite effect that passes the other filters. |
| `min_model_rsq` | `0.1` | Minimum finite target-model R² from `Pando::gof()`. |
| `require_padj` | `TRUE` | Require the coefficient table to contain valid adjusted P values. |

The Stage 1 evidence filter defines the Stage 3 supported metabolic-gene set. Positive and negative coefficients both count as regulatory evidence because the question is whether a gene has significant epigenetic regulatory support, not whether the inferred effect is activating or repressing.

## Metacell geometry and reproducibility

The four reduction fields determine the cell-level geometry used by SuperCell2:

| Parameter | Default | Meaning |
|---|---:|---|
| `rna_reduction` | `"pca"` | RNA reduction stored in `A@reductions`. A precomputed `"harmony"` reduction may be used instead. |
| `rna_dims` | `1:30` | RNA dimensions passed to SuperCell2. |
| `atac_reduction` | `"lsi"` | ATAC reduction stored in `A@reductions`. |
| `atac_dims` | `2:30` | ATAC dimensions passed to SuperCell2; LSI dimension 1 is excluded by default. |
| `seed` | `12345L` | Base random seed. For ordered condition strata, the internal seed is `seed + stratum_index - 1`. |
| `gamma` | `30L` | Approximate cells-per-metacell compression target. |
| `min_cells_per_stratum` | `100L` | Minimum cells required for a condition stratum. |
| `min_metacell_size` | `20L` | Metacells below this size are marked low-power. |
| `min_metacells_per_stratum` | `2L` | Minimum accepted number of metacells per condition stratum. |
| `overwrite` | `FALSE` | Set to `TRUE` when rebuilding checkpoints after changing cells, assays, reductions, dimensions, seed, gamma, or thresholds. |

Before running, verify that the requested reductions and dimensions exist:

```r
stopifnot(
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions),
  ncol(SeuratObject::Embeddings(A[["pca"]])) >= 30,
  ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
)
```

Using Harmony changes only the RNA neighbourhood geometry used to assign metacell membership; RegCompass still aggregates the original RNA and ATAC assay counts. Do not use a Harmony embedding that removed the biological condition effect being analysed.

Stage 2 records the exact reduction names, dimensions, embedding fingerprints, seed, gamma, and metacell thresholds in:

```r
result$metacell_data$pooled$cache_contract$analysis_args
```

## Fixed Stage 3 expansion

Stage 3 always performs exactly one ordered expansion pass:

```text
core subsystem
→ KEGG/Reactome reaction equivalence
→ master-Rhea reaction equivalence
```

## GPR-AND aggregation

`layer1_args$gpr_and_method` controls genes joined by a GPR AND relationship. Allowed values are `"min"`, `"median"`, and `"mean"`; the default is `"min"`. Isozyme OR branches are summed in the canonical Layer 1 calculation.

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
