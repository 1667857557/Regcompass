# Tutorial Level 1: minimal one-shot run

This tutorial uses RegCompassR 2.0.0 with the Pando 1.4.0
`ConditionGRNFit` contract. It assumes a paired-cell RNA+ATAC Seurat object.

## 1. Required object state

The object must contain:

- normalized RNA and ATAC assays;
- the condition and cell-type metadata columns used below;
- RNA PCA and ATAC LSI reductions when the default metacell geometry is used;
- paired RNA and ATAC measurements for the same cells.

```r
stopifnot(
  all(c("Group", "cell_type") %in% colnames(A@meta.data)),
  "pca" %in% names(A@reductions),
  "lsi" %in% names(A@reductions),
  ncol(SeuratObject::Embeddings(A[["pca"]])) >= 30,
  ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
)
```

RegCompass intersects GEM GPR genes with RNA-assay row names and passes that
complete metabolic target set to Pando. Do not override `genes` inside
`pando_infer_args`.

## 2. Prepare the GEM and medium

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
  scenario = "high_glucose",
  species = "human"
)
```

When `pfm` is omitted, RegCompass loads Pando's bundled `motifs` object. Human
runs also use the bundled hg38 phastCons plus SCREEN ccRE union unless
`pando_initiate_args$regions` is supplied.

## 3. Run the workflow

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
      method = "shared_baseline_condition_sparse",
      candidate_screen = "motif_domain",
      tf_cor = 0.1,
      peak_cor = 0,
      alpha = 0.5,
      condition_mix = 0.5,
      condition_weight = "equal",
      reference_condition = "Control",
      nlambda = 50L,
      nfolds = 5L,
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
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  layer1_args = list(
    regulatory_alpha = 1,
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
  layer2_workers = 30L,
  progress = TRUE
)
```

The one-shot workflow controls Pando execution through `upstream_workers`.
Do not add `parallel` or `BPPARAM` to `pando_infer_args`.

## 4. Why `motif_domain` is the default

Pando fits the interaction predictor

```text
TF_RNA × peak_ATAC
```

A useful interaction can coexist with weak marginal TF-target and peak-target
correlations. `candidate_screen = "motif_domain"` therefore retains the
motif/domain-supported dictionary and leaves edge selection to elastic-net
regularization. `pooled_within_condition` remains an optional marginal-screen sensitivity mode; its response-dependent screen makes its OOF score ineligible for confirmatory Layer 1 reliability (`q = 0`).

The directly comparable coefficient contract requires:

```r
method = "shared_baseline_condition_sparse"
condition_mix = 0.5
condition_weight = "equal"
scale = TRUE
```

Set `reference_condition` explicitly. Pando otherwise uses the first condition
level within each cell type.

## 5. Comparison support

For condition `c` and reference `r`, Pando exports:

```text
comparison_mask[e, c] = eligibility_mask[e, c] && eligibility_mask[e, r]
```

RegCompass uses `beta[e, c] - beta[e, r]` only where this mask is true. This
prevents a coefficient fixed to zero because an edge was non-estimable from
being interpreted as a biological loss.

```r
head(result$grn$tf_peak_gene_condition_effect_all[, c(
  "edge_id",
  "Group",
  "cell_type",
  "condition_estimate",
  "reference_estimate",
  "condition_effect",
  "comparable_to_reference"
)])

table(
  result$grn$tf_peak_gene_condition_effect_all$comparable_to_reference,
  useNA = "ifany"
)
```

## 6. Inspect the main outputs

```r
result$grn$condition_grn_fits
result$grn$condition_fit_status
result$grn$tf_peak_gene_condition
result$grn$tf_peak_gene_condition_effect
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$reaction_ranking
result$condition_contrast
result$merged_grn_meta_modules$merged_core_reactions
result$microcompass$model_cache_summary
```

Stage 1 provenance includes the actual candidate policy and comparison rule:

```r
result$grn$normalization_policy[c(
  "pando_candidate_screen",
  "comparison_support",
  "reference_condition",
  "coefficient_scale",
  "min_model_rsq"
)]
```

## 7. Mouse input

The bundled Pando regulatory regions are hg38 and must not be used for mouse
ATAC coordinates. Supply a build-matched `GRanges` object:

```r
library(BSgenome.Mmusculus.UCSC.mm10)

mouse_regions <- readRDS("mm10_regulatory_regions.rds")
stopifnot(methods::is(mouse_regions, "GenomicRanges"))

mouse_result <- rc_run_regcompass_one_shot(
  object = A_mouse,
  outdir = "RegCompass_mouse",
  genome = BSgenome.Mmusculus.UCSC.mm10,
  species = "mouse",
  condition_col = "Group",
  celltype_col = "cell_type",
  pando_args = list(
    pando_initiate_args = list(regions = mouse_regions),
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      reference_condition = "Control"
    )
  )
)
```

The region build must match both the ATAC peak coordinates and the genome used
for motif scanning.

## 8. Next step

Use [Tutorial 2](tutorial-02-stepwise-audit.md) when each stage should be saved,
validated, and restarted independently. API index: [functions.md](functions.md).


Use `sample_col = "sample_id"` for biological-sample provenance.

Set `cv_block_col = "sample_id"` for sample-blocked OOF validation.

SuperCell hard strata are exactly condition × broad cell type.
