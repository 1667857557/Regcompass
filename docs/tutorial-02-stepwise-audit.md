# Tutorial Level 2: stepwise run

Use this tutorial when each RegCompass stage should be run and saved independently.

## Configure cross-platform parallel backends

The explicit stepwise workflow accepts `BiocParallelParam` objects through `BPPARAM`. Define the upstream and Layer 2 backends once before Stage 1:

```r
library(BiocParallel)

upstream_workers <- 6L
layer2_workers <- 30L

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(
    workers = upstream_workers,
    type = "SOCK",
    progressbar = TRUE
  )
} else {
  MulticoreParam(
    workers = upstream_workers,
    progressbar = TRUE
  )
}

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(
    workers = layer2_workers,
    type = "SOCK",
    progressbar = TRUE
  )
} else {
  MulticoreParam(
    workers = layer2_workers,
    progressbar = TRUE
  )
}
```

Windows uses socket workers because fork-based `MulticoreParam` is unavailable. Linux and macOS use forked `MulticoreParam` workers. Reuse `upstream_bp` for Stage 1 and Stage 4, and use `layer2_bp` for Stage 5. The values `6L` and `30L` are examples rather than universal defaults; do not request more workers than the CPU and memory allocation available to the R process or batch job.

The one-shot runner does not require these objects. It accepts `upstream_workers` and `layer2_workers` directly and resolves the operating-system backend automatically.

## Stage 1: infer condition-comparable Pando evidence

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 100,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "shared_design_independent",
      candidate_screen = "condition_union",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 1,
      condition_weight = "equal",
      reference_condition = "Control",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
      parallel = FALSE
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)

step1$grn_result$sample_status
step1$grn_result$target_metabolic_genes
```

The Stage 1 runner arguments are ordered as shared inputs → motif/region policy → metadata/assays → Pando settings → execution controls.

When `pfm` is omitted, RegCompass loads `data("motifs", package = "Pando")` and passes `motifs` to `Pando::find_motifs()`. Supply `pfm = custom_motifs` only to override this default.

The candidate targets are all GEM GPR genes present in the RNA assay. Unless `pando_args$pando_initiate_args$regions` is supplied, the default regions are equivalent to:

```r
# Load the Pando data objects.
data("phastConsElements20Mammals.UCSC.hg38", package = "Pando")
data("SCREEN.ccRE.UCSC.hg38", package = "Pando")

# Human default.
human_regions <- union(
  phastConsElements20Mammals.UCSC.hg38,
  SCREEN.ccRE.UCSC.hg38
)

# Mouse default.
mouse_regions <- phastConsElements20Mammals.UCSC.hg38
```

Human uses phastCons plus SCREEN ccRE; mouse uses only `phastConsElements20Mammals.UCSC.hg38`. An explicit region object overrides either default.

### Stage 1 filter meanings

| Parameter | Default | Effect |
|---|---:|---|
| `min_cells` | `20L` | Minimum cells in every condition of a cell-type fit. |
| `method` | `"shared_design_independent"` | Independent condition coefficients in one comparable design. |
| `candidate_screen` | `"condition_union"` | Union only complete edges retained within a condition. |
| `condition_mix` | `1` | No cross-condition group penalty. |
| `condition_weight` | `"equal"` | Equal-condition loss convention. |
| `reference_condition` | first level | Baseline for explicit coefficient differences; set it explicitly for an audit. |
| `scale` | `TRUE` | Pooled scaling of each final TF-RNA × peak-ATAC edge predictor. |
| `min_abs_estimate` | `0` | Minimum absolute condition coefficient/reference contrast. |
| `min_model_rsq` | `0.1` | Minimum finite condition target-model R². |

The retained gene set is based on active condition coefficients. Coefficient
direction does not affect inclusion. The penalty path uses
`β_condition - β_reference`, the stored Pando edge transform, and the actual
metacell TF RNA × peak ATAC predictor. The regularized solver has no adjusted
P-value filter.

Audit the lossless contract and transform exports before Stage 2:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition_effect

read.delim(gzfile(file.path(
  "RegCompass_steps/01_grn",
  "pando_edge_predictor_transforms.tsv.gz"
)))
```

## Stage 2: construct condition-level metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  sample_col = NULL,
  condition_col = "Group",
  celltype_col = "cell_type",
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
  )
)

step2$pooled$metacell_meta
step2$pooled$cache_contract$analysis_args
```

The Stage 2 defaults are:

```r
rna_reduction = "pca"
rna_dims = 1:30
atac_reduction = "lsi"
atac_dims = 2:30
gamma = 30L
seed = 12345L
min_cells_per_stratum = 100L
min_metacell_size = 20L
min_metacells_per_stratum = 2L
```

`rna_reduction` and `atac_reduction` must name reductions already present in `A@reductions`; the selected dimensions must exist. LSI dimension 1 is excluded by default. The base seed is incremented deterministically by condition-stratum order:

```text
internal_seed = seed + stratum_index - 1
```

Before running the default geometry, verify:

```r
stopifnot(
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions),
  ncol(SeuratObject::Embeddings(A[["pca"]])) >= 30,
  ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
)
```

A precomputed Harmony embedding may replace PCA:

```r
metacell_args = list(
  rna_reduction = "harmony",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30,
  seed = 12345L,
  overwrite = TRUE
)
```

Harmony affects only the RNA neighbourhood geometry used for membership construction; original assay counts are still aggregated. Do not use a Harmony embedding that removed the biological condition contrast. Changing cells, assay matrices, reduction names, dimensions, embedding values, seed, gamma, or metacell thresholds invalidates Stage 2 checkpoints; use `overwrite = TRUE` to rebuild.

## Stage 3: construct complete-GPR biological meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules"
)
```

For each `condition × cell type`, Stage 3 performs the following operations exactly once:

```text
active Pando TF–peak–GEM-target rows
→ unique supported metabolic target genes
→ complete-GPR core reactions
→ all reactions in core-reaction subsystems
→ direct KEGG/Reactome reaction-equivalence expansion
→ direct master-Rhea reaction-equivalence expansion
→ biological meta-module
```

A positive or negative Pando coefficient both count as regulatory evidence. A reaction is core only when at least one complete GPR AND branch is contained in the supported target-gene set. Partial complexes remain diagnostic and do not anchor expansion.

A custom subsystem mapping remains possible through:

```r
step3_custom <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules_custom",
  meta_module_args = list(
    subsystem_table = custom_subsystem_table
  )
)
```

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
table(step3$condition_modules$reaction_membership$inclusion_stage)
step3$condition_modules$meta_module_summary$expansion_policy

catalogue <- step3$merged_modules
catalogue$merged_core_reactions
catalogue$merged_reaction_membership
```

## Stage 4: calculate integrated RNA+ATAC reaction support

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  gpr_and_method = "min",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

`gpr_and_method` accepts the COMPASS functions `"min"`, `"median"`, and `"mean"`. RegCompass defaults to `"min"`, so the least-supported required subunit limits a multi-gene GPR complex. The canonical isozyme OR rule remains additive.

The selected rule is recorded in:

```r
step4$capacity_params$and_method
step4$evidence_formula
```

## Stage 5: build the medium-constrained model and score reactions

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2",
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
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

Stage 5 first applies the selected medium and performs global FASTCORE completion to construct the union GEM. `completion_time_limit` applies only to this construction phase. The completed union GEM is then cached and reused for directional scoring; scoring LPs have no time-limit parameter. See [medium presets](medium-presets.md) for available presets and custom media.

```r
step5$model_cache_summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "n_reactions",
  "file"
)]
```

## Stage 6: assemble annotated results

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results",
  species = "human"
)
```

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
```
