# Tutorial Level 2: stepwise run with audit gates

Use the stepwise workflow when intermediate objects must be inspected or restarted. Each stage writes a classed RDS object. RegCompassR 1.8.2 verifies workflow parameters, GEM fingerprints, and ordered metacell IDs before accepting downstream inputs.

## Setup

```r
library(RegCompassR)
library(Pando)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")
condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"

upstream_bp <- MulticoreParam(workers = 16L, progressbar = TRUE)
layer2_bp <- MulticoreParam(workers = 12L, progressbar = TRUE)

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "physiologic",
  species = "human"
)
```

Before starting R on Linux:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## Stage 1: condition-by-cell-type GRNs

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn",
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
  BPPARAM = upstream_bp
)

status1 <- step1$grn_result$sample_status
stopifnot(
  inherits(step1, "regcompass_grn_step"),
  all(status1$status == "ok"),
  all(status1$n_significant_edges > 0),
  nzchar(step1$gem_fingerprint)
)
```

RNA is normalized once. ATAC TF-IDF uses all conditions within each cell type as the reference. Keep Pando's inner `parallel = FALSE` when the outer stage is parallel.

## Stage 2: condition-level metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = condition_col,
  celltype_col = celltype_col,
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  )
)

meta2 <- step2$pooled$metacell_meta
stopifnot(
  inherits(step2, "regcompass_metacell_step"),
  !any(meta2$dominant_celltype_tied %in% TRUE),
  setequal(colnames(step2$metacell_object), meta2$metacell_id)
)
```

Condition is the only hard stratum. `celltype_col` guides SuperCell2 before aggregation and is audited again from member-cell composition. Stage 2 does not accept workflow-level `BPPARAM`.

## Stage 3: core reactions and meta-modules

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
      workers = 16L,
      backend = "multicore"
    )
  )
)

stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  all(step3$group_coverage$coverage_complete),
  nrow(step3$global_modules$global_core_reactions) > 0,
  nrow(step3$global_modules$global_reaction_membership) > 0,
  identical(step3$gem_fingerprint, step1$gem_fingerprint)
)
```

Biological modules contain complete-GPR cores, same-subsystem reactions, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Local FASTCORE adds only reactions needed for feasibility.

## Stage 4: integrated reaction expression

```r
step4 <- rc_regcompass_step_layer1(
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_steps/04_layer1",
  regulatory_alpha = 1,
  tau = 0.20,
  parallel = TRUE,
  BPPARAM = upstream_bp
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
```

## Stage 5: original core LP scoring

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
  BPPARAM = layer2_bp
)

stopifnot(
  inherits(step5, "regcompass_layer2_step"),
  any(step5$evaluated),
  identical(colnames(step5$penalty), colnames(step4$reaction_expression)),
  identical(step5$workflow_params, step4$workflow_params),
  identical(step5$gem_fingerprint, step4$gem_fingerprint),
  all(file.exists(step5$model_cache_summary$file))
)
```

`penalty` is the primary output. `score` is a within-target relative rank. The union-GEM cache must be retained for restart and expanded-target analysis.

## Optional Stage 5b: score related reactions in the same union GEM

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/05b_expanded_targets",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)

stopifnot(
  inherits(expanded, "regcompass_target_union_step"),
  all(expanded$expanded_scoring_targets$score_target),
  expanded$microcompass$params$structural_model_reused_exactly,
  all(expanded$microcompass$model_cache_summary$reused_without_rebuilding)
)
```

This second pass changes only the LP target. It does not rebuild FASTCORE support, stoichiometry, medium bounds, or the union reaction set.

## Stage 6: assemble results

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/06_results"
)

stopifnot(
  identical(result$version, "1.8.2"),
  identical(result$schema_version, "regcompass_grn_first_v2"),
  nrow(result$reaction_ranking) > 0,
  nrow(result$reaction_catalog) > 0,
  nrow(result$reaction_evidence) > 0,
  identical(result$gem_fingerprint, step5$gem_fingerprint)
)
```

## Stop conditions

Do not continue when a required GRN group fails, a metacell has a tied dominant label, GRN/metacell coverage is incomplete, no complete-GPR core remains, stage classes or fingerprints differ, unit order changes, the solver is unavailable, or no LP target is evaluated.
