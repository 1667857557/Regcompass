# Tutorial Level 2: audit the saved stage workflow

Use the canonical complete workflow for computation and the saved classed stage objects for inspection, restart planning, and audit gates. RegCompassR 1.8.3 automatically applies the two worker layers, selects the platform backend, and releases each parallel stage pool immediately after completion or failure.

## Setup

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")
condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"

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

Before starting R on Linux:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## Run all stages with automatic worker management

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
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
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  layer1_args = list(
    top_k_neighbors = 5,
    min_shared_tfs = 1,
    min_tf_jaccard = 0,
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      strict = TRUE,
      time_limit = 300
    )
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  upstream_workers = 6L,
  layer2_workers = 30L,
  progress = TRUE
)
```

The complete run writes every intermediate wrapper RDS. The canonical interface does not expose `parallel_backend`; Windows uses SOCK workers, Linux/macOS use multicore workers, and either worker count set to one resolves that layer to serial execution.

## Load the saved stage objects

```r
step1 <- readRDS(
  "RegCompass_result/01_single_cell_grn/step_grn.rds"
)
step2 <- readRDS(
  "RegCompass_result/02_condition_metacells/step_metacells.rds"
)
step3 <- readRDS(
  "RegCompass_result/03_meta_modules/step_meta_modules.rds"
)
step4 <- readRDS(
  "RegCompass_result/04_layer1/step_layer1.rds"
)
step5 <- readRDS(
  "RegCompass_result/05_layer2/step_layer2.rds"
)
result <- readRDS(
  "RegCompass_result/06_results/regcompass_result.rds"
)
```

Use the wrapper RDS files above rather than compact inspection artifacts.

## Audit Stage 1: condition-by-cell-type GRNs

```r
status1 <- step1$grn_result$sample_status

stopifnot(
  inherits(step1, "regcompass_grn_step"),
  all(status1$status == "ok"),
  all(status1$n_significant_edges > 0),
  nzchar(step1$gem_fingerprint)
)

status1
```

RNA is normalized once. ATAC TF-IDF uses all conditions within each cell type as the reference. Pando's inner `parallel = FALSE` is retained because condition-by-cell-type groups are distributed by the upstream worker layer.

## Audit Stage 2: condition-level metacells

```r
meta2 <- step2$pooled$metacell_meta

stopifnot(
  inherits(step2, "regcompass_metacell_step"),
  !any(meta2$dominant_celltype_tied %in% TRUE),
  setequal(
    colnames(step2$metacell_object),
    meta2$metacell_id
  )
)

table(
  meta2[[condition_col]],
  meta2[[celltype_col]]
)
```

Condition is the only hard stratum. `celltype_col` guides SuperCell2 before aggregation and is audited again from member-cell composition. Stage 2 has no workflow-level parallel worker pool.

## Audit Stage 3: core reactions and meta-modules

```r
stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  all(step3$group_coverage$coverage_complete),
  nrow(step3$global_modules$global_core_reactions) > 0,
  nrow(step3$global_modules$global_reaction_membership) > 0,
  identical(step3$gem_fingerprint, step1$gem_fingerprint)
)

step3$condition_modules$local_fastcore_summary
```

Biological meta-modules contain complete-GPR cores, same-subsystem reactions, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Local FASTCORE adds only reactions required for feasibility. Its workers are assigned from `upstream_workers` and released before Layer 1 begins.

## Audit Stage 4: integrated reaction expression

```r
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

Stage 4 receives a fresh upstream worker pool. That pool is stopped and garbage-collected immediately after the stage returns.

## Audit Stage 5: original core LP scoring

```r
stopifnot(
  inherits(step5, "regcompass_layer2_step"),
  any(step5$evaluated),
  identical(
    colnames(step5$penalty),
    colnames(step4$reaction_expression)
  ),
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

`penalty` is the primary output. `score` is a within-target relative rank. Layer 2 uses its own `layer2_workers` pool; no upstream pool remains active at this point.

## Optional direct database-linked target scoring

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_result/05b_expanded_targets",
  core_genes = c("GCLC", "GCLM", "GSS", "GSR", "G6PD", "PGD"),
  gene_match = "complete_gpr",
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  progress = TRUE
)
```

The selected cores are not scored again. The second pass includes only non-core reactions directly sharing KEGG, Reactome, or master-Rhea identifiers with the selected cores.

## Audit Stage 6 and execution provenance

```r
stopifnot(
  identical(result$version, "1.8.3"),
  identical(result$schema_version, "regcompass_grn_first_v2"),
  nrow(result$reaction_ranking) > 0,
  nrow(result$reaction_catalog) > 0,
  nrow(result$reaction_evidence) > 0,
  identical(result$gem_fingerprint, step5$gem_fingerprint),
  identical(result$params$upstream_workers, 6L),
  identical(result$params$layer2_workers, 30L),
  identical(
    result$params$parallel_worker_lifecycle,
    "stage_scoped_create_start_stop_release_full_gc"
  )
)

result$params$parallel_backend_resolved
result$params$parallel_stage_groups
result$timing$stages
```

## Stop conditions

Do not continue when a required GRN group fails, a metacell has a tied dominant label, GRN/metacell coverage is incomplete, no complete-GPR core remains, stage classes or fingerprints differ, unit order changes, the solver is unavailable, or no LP target is evaluated.
