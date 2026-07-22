# Tutorial Level 2: stepwise run with audit gates

Use this level when every stage must be inspected before the next stage starts. It follows the same analysis design as Level 1 but exposes the intermediate objects and files.

Complete the installation and input checks in [Level 1](tutorial-01-quick-start.md) first. Use [Level 3](tutorial-03-advanced-restart.md) for restart, alternative media, solver selection, and detailed resource controls.

## Stage map

| Stage | Primary input | Primary output | Parallel unit on Linux | Required gate |
|---|---|---|---|---|
| 1. GRN | single-cell Seurat object, GEM, motifs, genome | condition × cell-type Pando GRNs | one condition × cell-type group per worker | every required group is `ok` and has significant edges |
| 2. Metacells | original Seurat object and complete annotation label | label-guided metacells within condition plus composition audit | not controlled by the workflow `BPPARAM` | no ambiguous dominant-cell-type ties; inspect purity and mixing |
| 3. Meta-modules | Stages 1-2, GEM | core reactions and expanded modules | one local FASTCORE meta-module completion per worker | GRN/metacell coverage is complete; core reactions exist |
| 4. Layer 1 | metacells, modules, GEM | reaction-expression matrix | GPR/reaction-capacity calculations | columns align exactly to metacell metadata |
| 5. Layer 2 | Layer 1, global module, medium | directional LP scores | one shared-model × metacell task per worker | targets were evaluated and feasible targets exist |
| 6. Results | Stages 1-5 | rankings and condition contrasts | serial assembly | outputs retain condition-specific and global modules |

## Common Linux multicore setup

This tutorial assumes a normal Linux host, not Windows. Use explicit `MulticoreParam` objects so worker allocation is visible and reproducible.

```r
library(RegCompassR)
library(Pando)
library(Signac)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")

condition_col <- "dataset"
celltype_col <- "epithelial_or_stem"
metacell_label_col <- celltype_col

upstream_workers <- 16L
layer2_workers <- 12L

upstream_bp <- BiocParallel::MulticoreParam(
  workers = upstream_workers,
  progressbar = TRUE
)

layer2_bp <- BiocParallel::MulticoreParam(
  workers = layer2_workers,
  progressbar = TRUE
)

gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "high_glucose",
  species = "human"
)
rc_validate_gem(gem)
```

The same paired-cell Seurat object `A` is passed independently to Stages 1 and 2. Do not pass the internally normalized Stage 1 object into Stage 2.

Set numerical-library thread counts to one before starting R when the server uses OpenBLAS or MKL. This prevents every forked worker from starting additional BLAS threads:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

## Stage 1: condition × cell-type single-cell GRNs

### Input

- original paired-cell RNA+ATAC object `A`;
- validated `gem`;
- `motifs`, not `motif2tf`;
- genome matching ATAC peak coordinates;
- complete condition and cell-type metadata.

### Run with group-level multicore parallelism

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
```

RNA is normalized across all cells. ATAC TF-IDF uses all conditions within each cell type as the shared reference. RegCompass distributes the independent condition × cell-type subsets across the outer workers.

Keep `pando_infer_args$parallel = FALSE`. Activating both the outer RegCompass workers and Pando's internal workers creates nested parallelism and usually oversubscribes the Linux host.

### Gate before Stage 2 or 3

```r
status1 <- step1$grn_result$sample_status
status1

stopifnot(
  all(status1$status == "ok"),
  all(status1$n_significant_edges > 0),
  nrow(step1$grn_result$tf_peak_gene_significant) > 0
)

step1$grn_result$pando_installation
head(step1$grn_result$tf_peak_gene_significant)
```

Do not continue when a required group failed, was skipped for too few cells, or has no significant edge.

### Files

- `single_cell_grn.rds`;
- `step_grn.rds`;
- `pando_group_status.tsv.gz`;
- `pando_tf_peak_gene_all.tsv.gz`;
- `pando_tf_peak_gene_significant.tsv.gz`;
- optional `pando_objects/*.rds`.

## Stage 2: label-aware, condition-only metacells

### Input

The original `A`. Condition is the only hard metacell stratum. The existing
cell-type annotation is supplied to SuperCell2 as a construction label so cells
of different annotated types are not merged indiscriminately. Cell type and
sample remain excluded from the hard stratum definition.

### Run

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_metacells",
  condition_col = condition_col,
  celltype_col = celltype_col,
  label_col = metacell_label_col,
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 75,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  )
)
```

Stage 2 does not use the workflow `BPPARAM` shown above. Do not insert `BPPARAM = TRUE` into `metacell_args`. `label_col` must name a complete annotation column and is passed to SuperCell2 before aggregation. Each metacell receives a dominant member-cell type after construction. Purity, mixed-cell-type status, and the full composition are retained. An exact dominant-cell-type tie stops the workflow.

### Gate before Stage 3

```r
meta2 <- step2$pooled$metacell_meta

stopifnot(
  nrow(meta2) > 0,
  !anyNA(meta2[[condition_col]]),
  !anyNA(meta2[[celltype_col]]),
  !any(meta2$dominant_celltype_tied %in% TRUE),
  setequal(colnames(step2$metacell_object), meta2$metacell_id)
)

table(meta2[[condition_col]], meta2[[celltype_col]])
summary(meta2$dominant_celltype_fraction)
```

### Files

- `step_metacells.rds`;
- `merged_metacell_object.rds`;
- `metacell_metadata.tsv.gz`;
- `metacell_membership.tsv.gz`;
- `metacell_celltype_composition.tsv.gz`;
- `metacell_celltype_summary.tsv.gz`.

## Stage 3: core reactions and meta-modules

### Input

Stage 1, Stage 2, and the same GEM used for Stage 1.

### Run local FASTCORE completion in parallel

The projection, core-reaction mapping, and database expansion are deterministic setup operations. The expensive local FASTCORE completion is distributed by meta-module. The worker configuration is placed inside `local_fastcore_args` because that is the parallel substage.

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
      backend = "multicore"
    )
  )
)
```

The stage:

1. validates bidirectional condition × cell-type coverage between GRNs and metacells;
2. projects the Pando network onto metabolic genes;
3. maps complete-GPR core reactions;
4. adds reactions from core-reaction subsystems;
5. adds reactions sharing KEGG, Reactome, or identical master-Rhea identifiers;
6. completes separate local meta-modules with FASTCORE in parallel.

The parent GEM is prepared once before the worker loop. Each worker completes a different `sample_id × module_id`, so concurrent workers do not alter one another's biological module definition.

### Gate before Stage 4

```r
coverage3 <- step3$group_coverage
core3 <- step3$condition_modules$core_gene_reaction
membership3 <- step3$condition_modules$reaction_membership
fastcore3 <- step3$condition_modules$local_fastcore_summary

stopifnot(
  all(coverage3$coverage_complete),
  nrow(core3) > 0,
  nrow(membership3) > 0,
  all(core3$reaction_id %in% colnames(gem$S)),
  all(membership3$reaction_id %in% colnames(gem$S)),
  all(fastcore3$parallel_task == "local_fastcore_by_meta_module"),
  all(fastcore3$parallel_workers == upstream_workers)
)

coverage3
head(core3)
head(membership3)
unique(fastcore3[, c(
  "parallel_task", "parallel_backend", "parallel_workers"
)])
```

### Files

- `grn_metacell_group_coverage.tsv.gz`;
- `core_gene_reaction.tsv.gz`;
- `meta_module_reactions.tsv.gz`;
- `condition_meta_modules.rds`;
- `global_meta_modules.rds`;
- `local_fastcore/` diagnostics and optional per-module models.

## Stage 4: RNA+ATAC reaction expression

### Run with GPR/reaction parallelism

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
```

### Gate before Stage 5

```r
stopifnot(
  nrow(step4$reaction_expression) > 0,
  ncol(step4$reaction_expression) > 0,
  identical(
    colnames(step4$reaction_expression),
    as.character(step4$unit_meta$pool_id)
  ),
  all(is.finite(step4$reaction_expression))
)

dim(step4$reaction_expression)
head(step4$gpr_diagnostics)
step4$capacity_params[c("parallel", "bpparam_class")]
```

`step_layer1.rds` contains RNA support, the ATAC regulatory modifier, integrated gene support, GPR diagnostics, and reaction expression.

## Stage 5: directional COMPASS-like scoring

### Run LP tasks with a separate worker pool

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
```

Layer 2 first builds and caches each structural GEM. It then distributes the independent `shared model × metacell` LP tasks. Model-cache construction itself remains serial so workers do not race while writing the same cache file.

Use a smaller `layer2_workers` value than the CPU count when each GEM is large. Parallel LP tasks can become memory-bound before they become CPU-bound.

The selected solver is checked before medium-constrained model construction. A missing solver package is reported separately from true biological infeasibility.

### Gate before Stage 6

```r
stopifnot(
  any(step5$evaluated),
  nrow(step5$lp_diagnostics) > 0,
  nrow(step5$model_cache_summary) > 0,
  identical(step5$params$parallel_task, "shared_model_by_metacell")
)

table(step5$evaluated)
table(step5$feasible)
head(step5$lp_diagnostics)
step5$model_cache_summary
```

A run can contain biologically blocked target directions, but a result with no evaluated target or no feasible target requires inspection before interpretation.

### Files

- `step_layer2.rds`;
- exported score, penalty, `vmax`, feasibility, and diagnostic tables;
- persistent `model_cache/`.

## Stage 6: final result assembly

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
```

### Final gate

```r
stopifnot(
  nrow(result$reaction_ranking) > 0,
  identical(result$version, "1.8.1"),
  !is.null(result$condition_grn_meta_modules),
  !is.null(result$global_grn_meta_modules)
)

head(result$reaction_ranking)
head(result$condition_summary)
head(result$condition_contrast)
```

Final files are `step_comparison.rds` and `regcompass_result.rds`.

## Stop conditions

Stop and diagnose rather than continuing when any of these occur:

- a condition × cell-type Pando group is not `ok`;
- a metacell has an ambiguous dominant-cell-type tie;
- GRN/metacell group coverage is incomplete;
- no complete-GPR core reaction remains;
- Layer 1 reaction-expression columns do not align to metacell metadata;
- Layer 2 evaluated no target;
- solver installation is reported as missing.
