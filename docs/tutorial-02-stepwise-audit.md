# Tutorial Level 2: true stepwise run with audit gates

Use this tutorial when every RegCompass stage must be executed, inspected, saved, and restarted independently. Unlike the one-shot tutorial, this workflow calls each public stage function directly.

RegCompassR 1.8.3 validates stage classes, workflow metadata, GEM fingerprints, core-target provenance, and ordered metacell IDs before accepting downstream objects.

## Setup

```r
library(RegCompassR)
library(Pando)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"
upstream_workers <- 6L
layer2_workers <- 30L

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

Use `Pando::motifs` as `pfm`; `motif2tf` is not a PFM/PWM collection.

Before starting R on Linux, keep numerical-library threads at one because the outer BiocParallel layer distributes independent tasks:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Create stage-specific worker backends:

```r
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

For serial troubleshooting, set `parallel = FALSE`, omit `BPPARAM`, and set local FASTCORE `parallel = FALSE`.

## Stage 1: infer condition-by-cell-type GRNs

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_single_cell_grn",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = condition_col,
  celltype_col = celltype_col,
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
  parallel = TRUE,
  BPPARAM = upstream_bp,
  progress = TRUE
)

status1 <- step1$grn_result$sample_status
stopifnot(
  inherits(step1, "regcompass_grn_step"),
  all(status1$status == "ok"),
  all(status1$n_significant_edges > 0),
  nzchar(step1$gem_fingerprint)
)
```

RNA is normalized once. ATAC TF-IDF uses all conditions within each cell type as the reference, and Pando is fitted separately for each `condition × cell type` group. Pando's inner `parallel = FALSE` prevents nested parallelism.

## Stage 2: construct condition-level metacells

By default, SuperCell2 uses RNA `pca` dimensions 1:30 and ATAC `lsi` dimensions 2:30. Set the following switch to `TRUE` only when an existing Harmony reduction should define the RNA geometry used for metacell construction.

```r
use_harmony_for_metacells <- FALSE

metacell_embedding_args <- if (use_harmony_for_metacells) {
  stopifnot(
    "harmony" %in% names(A@reductions),
    "lsi" %in% names(A@reductions),
    ncol(SeuratObject::Embeddings(A[["harmony"]])) >= 30,
    ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
  )
  list(
    rna_reduction = "harmony",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30
  )
} else {
  list()
}

step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  sample_col = NULL,
  condition_col = condition_col,
  celltype_col = celltype_col,
  fragment_files = FALSE,
  metacell_args = c(
    metacell_embedding_args,
    list(
      gamma = 30,
      min_cells_per_stratum = 500,
      min_metacell_size = 10
    )
  ),
  progress = TRUE
)

meta2 <- step2$pooled$metacell_meta
stopifnot(
  inherits(step2, "regcompass_metacell_step"),
  identical(step2$params$metacell_args$gamma, 30L),
  identical(
    step2$pooled$cache_contract$analysis_args$rna_reduction,
    if (use_harmony_for_metacells) "harmony" else "pca"
  ),
  !any(meta2$dominant_celltype_tied %in% TRUE),
  setequal(colnames(step2$metacell_object), meta2$metacell_id)
)

table(meta2[[condition_col]], meta2[[celltype_col]])
```

Condition is the only hard stratum. `celltype_col` is passed to SuperCell2 as a construction label and is audited again from member-cell composition. `sample_col` is optional provenance and is not used for balancing.

Harmony affects only the RNA embedding used to define cell-to-metacell geometry. RNA and ATAC metacell counts remain sums from the original assays. Switching between PCA and Harmony, changing dimensions, or recomputing either reduction changes the Stage 2 cache contract; existing checkpoints must be rebuilt with `overwrite = TRUE`.

## Stage 3: build core reactions and meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      strict = TRUE,
      time_limit = 300,
      parallel = TRUE,
      workers = upstream_workers,
      backend = "auto"
    )
  ),
  progress = TRUE
)

stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  all(step3$group_coverage$coverage_complete),
  nrow(step3$global_modules$global_core_reactions) > 0,
  nrow(step3$global_modules$global_reaction_membership) > 0,
  identical(step3$gem_fingerprint, step1$gem_fingerprint)
)

step3$condition_modules$local_fastcore_summary
```

Biological meta-modules contain complete-GPR cores, same-subsystem reactions, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Local FASTCORE adds separately labelled reactions required for feasibility.

## Stage 4: calculate integrated RNA+ATAC reaction expression

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20,
  parallel = TRUE,
  BPPARAM = upstream_bp,
  progress = TRUE
)

stopifnot(
  inherits(step4, "regcompass_layer1_step"),
  identical(
    colnames(step4$reaction_expression),
    as.character(step4$unit_meta$pool_id)
  ),
  identical(step4$workflow_params, step3$workflow_params),
  identical(step4$gem_fingerprint, step3$gem_fingerprint)
)

dim(step4$reaction_expression)
head(step4$gpr_diagnostics)
```

## Stage 5: run directional COMPASS-like LP scoring

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
    time_limit = 60
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp,
  progress = TRUE
)

stopifnot(
  inherits(step5, "regcompass_layer2_step"),
  any(step5$evaluated),
  identical(colnames(step5$penalty), colnames(step4$reaction_expression)),
  identical(step5$workflow_params, step4$workflow_params),
  identical(step5$gem_fingerprint, step4$gem_fingerprint),
  all(file.exists(step5$model_cache_summary$file))
)

table(
  evaluated = as.vector(step5$evaluated),
  feasible = as.vector(step5$feasible),
  useNA = "ifany"
)
```

`penalty` is the primary output. Lower values indicate stronger compatibility between the required target flux and integrated evidence. `score` is a within-target relative rank, not a probability or measured flux.

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
  species = "human",
  progress = TRUE
)

stopifnot(
  identical(result$version, "1.8.3"),
  identical(result$schema_version, "regcompass_grn_first_v2"),
  nrow(result$reaction_ranking) > 0,
  nrow(result$reaction_catalog) > 0,
  nrow(result$reaction_evidence) > 0,
  identical(result$gem_fingerprint, step5$gem_fingerprint),
  identical(result$params$execution_mode, "stepwise")
)
```

Each stage writes its own classed RDS and `step_timing.tsv`, so a downstream stage can be restarted without recomputing valid upstream objects.

Stop the user-created workers after all stages finish:

```r
bpstop(upstream_bp)
bpstop(layer2_bp)
gc(full = TRUE)
```

## Optional Stage 5b: select anchors and remap linked reactions

A targeted second pass may be run from the completed Stage 3–5 objects. The selected genes or reaction IDs must resolve to previous global core reactions.

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/05b_targeted_reactions",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp,
  progress = TRUE
)

stopifnot(
  inherits(expanded, "regcompass_target_union_step"),
  expanded$microcompass$params$structural_model_reused_exactly,
  identical(
    expanded$microcompass$params$target_scope,
    "direct_kegg_reactome_master_rhea_noncore_only"
  )
)
```

The complete gene- and reaction-selection tutorial is [Tutorial Level 4](tutorial-04-targeted-reaction-remapping.md).

## Compare conditions after Stage 6

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = condition_col,
  celltype_col = celltype_col,
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  p_adjust_scope = "celltype_contrast_medium"
)

condition_stats$omnibus
condition_stats$pairwise
```

The complete statistics and plotting tutorial is [Tutorial Level 5](tutorial-05-condition-differential-analysis.md).

## Stop conditions

Do not continue when a required GRN group fails, a metacell has a tied dominant label, GRN/metacell coverage is incomplete, no complete-GPR core remains, stage classes or GEM fingerprints differ, scoring-unit order changes, the LP solver is unavailable, or no Layer 2 target is evaluated.
