# Tutorial Level 3: restart, sensitivity runs, Linux parallelism, and diagnostics

Use this level after completing the [Level 2 audit workflow](tutorial-02-stepwise-audit.md). It explains which stages must be rerun after a parameter change, how to allocate Linux workers, and how to separate installation, structural, and biological failures.

## 1. Restart from saved stage objects

For a stepwise run:

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
result <- readRDS("RegCompass_steps/06_results/regcompass_result.rds")
```

For a one-shot run, the corresponding stage directories are:

```text
01_single_cell_grn
02_condition_metacells
03_meta_modules
04_layer1
05_layer2
06_results
```

The one-shot workflow writes the final result both to
`RegCompass_result/regcompass_result.rds` and to
`RegCompass_result/06_results/regcompass_result.rds`. The two files contain the
same one-shot execution metadata, so either location is safe for a restart.

Always load the stage wrapper, such as `step_grn.rds` or `step_metacells.rds`, when a downstream function requires the class and workflow parameters. The compact files such as `single_cell_grn.rds` are useful for inspection but are not substitutes for the stage wrapper.

## 2. Minimal rerun matrix

| Changed item | Earliest stage to rerun | Reusable upstream stages |
|---|---:|---|
| Pando `tf_cor`, `peak_cor`, model method, minimum cells | 1 | none for GRN-derived modules; Stage 2 may be reused if unchanged |
| metacell `gamma`, minimum stratum size, minimum metacell size | 2 | Stage 1 |
| SuperCell2 `metacell_label_col` / Stage 2 `label_col` | 2 | Stage 1 |
| GRN projection or meta-module expansion settings | 3 | Stages 1-2 |
| local FASTCORE solver, strictness, or support limits | 3 | Stages 1-2 |
| `regulatory_alpha`, GPR `tau`, RNA half-saturation | 4 | Stages 1-3 |
| medium scenario, solver, `omega`, target direction, `model_mode` | 5 | Stages 1-4 |
| worker count or parallel backend only | same numerical stage | all completed upstream stages |
| display or final assembly only | 6 | Stages 1-5 |

When Stage 1 or Stage 2 changes, rerun Stage 3 because the bidirectional GRN/metacell group-coverage contract must be revalidated.

Changing only worker count or backend should not change the mathematical result. It does require rerunning the affected stage if you want the new resource setting to take effect.

## 3. Linux process and thread controls

RegCompass uses BiocParallel for independent R tasks. On a normal Linux host, `MulticoreParam` uses forked worker processes.

Before launching R, limit numerical-library threads so each forked worker does not start additional OpenMP, OpenBLAS, or MKL threads:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

An optional default worker count can be set for calls that use automatic discovery:

```bash
export REGCOMPASS_WORKERS=16
```

Explicit R objects are preferred in an analysis script:

```r
library(BiocParallel)

upstream_workers <- 16L
layer2_workers <- 12L

upstream_bp <- MulticoreParam(
  workers = upstream_workers,
  progressbar = TRUE
)

layer2_bp <- MulticoreParam(
  workers = layer2_workers,
  progressbar = TRUE
)
```

Use fewer workers than physical memory permits. Layer 2 often becomes memory-limited because several solver tasks can hold a GEM and penalty vectors at the same time.

## 4. Parallel units by workflow stage

| Stage | Parallel unit | Control | Important limitation |
|---|---|---|---|
| 1. Pando GRN | condition × cell-type group | `parallel`, `BPPARAM` | keep Pando inner `parallel = FALSE` |
| 2. Metacells | no workflow-level BiocParallel unit | SuperCell/metacell arguments only | do not pass `BPPARAM = TRUE` |
| 3. Meta-modules | local FASTCORE completion per `sample_id × module_id` | `local_fastcore_args$parallel`, `$workers`, `$backend`, or `$BPPARAM` | projection and identifier expansion are serial setup |
| 4. Layer 1 | GPR/reaction-capacity work | `parallel`, `BPPARAM` | ATAC modifier setup is not separately forked |
| 5. Layer 2 | unique shared-model × metacell task | `parallel`, `BPPARAM` | structural model-cache construction remains serial |
| 6. Results | serial assembly | none | no benefit from extra workers |

The workflow deliberately avoids automatic nested parallelism. A condition × cell-type group is one outer Stage 1 task, so `pando_infer_args$parallel` should remain `FALSE`.

## 5. One-shot Linux multicore run

`upstream_workers` is used for Stage 1 group-level Pando, Stage 3 local FASTCORE completion, and Stage 4 reaction capacity. `layer2_workers` is reserved for the larger number of Layer 2 LP tasks.

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result_linux",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  gem = gem,
  medium_scenarios = medium_scenarios,
  condition_col = condition_col,
  celltype_col = celltype_col,
  metacell_label_col = celltype_col,
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
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      strict = TRUE,
      parallel = TRUE
    )
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  upstream_workers = 16L,
  layer2_workers = 12L,
  parallel_backend = "multicore"
)
```

`parallel_backend = "auto"` also selects multicore on an ordinary non-container Linux host. Explicit `"multicore"` is preferable in a reproducible Linux script. In Docker or another detected container, `"auto"` selects a socket-based SnowParam backend instead.

`metacell_label_col` is evaluated in Stage 2. Changing it therefore requires
rerunning Stage 2 and every downstream stage, while an unchanged Stage 1 GRN may
be reused. The selected label column must be complete. It guides SuperCell2
before aggregation but does not replace the Stage 2 membership and purity audit.

## 6. Stepwise Linux multicore run

### Stage 1

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

### Stage 3 local FASTCORE

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules",
  layer1_args = list(
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

A preconstructed parameter object is also accepted at the internal FASTCORE boundary:

```r
layer1_args = list(
  local_fastcore = TRUE,
  local_fastcore_args = list(
    solver = "highs",
    parallel = TRUE,
    BPPARAM = upstream_bp
  )
)
```

For restartable stage objects, `workers` plus `backend` is usually cleaner than storing a live BiocParallel object inside the argument list.

### Stage 4

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

### Stage 5

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

## 7. Verify that parallel settings were applied

Stage 3 writes its worker policy into the local FASTCORE summary:

```r
unique(step3$condition_modules$local_fastcore_summary[, c(
  "parallel_task",
  "parallel_backend",
  "parallel_workers"
)])
```

Expected Linux multicore values include:

```text
parallel_task     = local_fastcore_by_meta_module
parallel_backend  = MulticoreParam
parallel_workers  = 16
```

Stage 4 records its backend class:

```r
step4$capacity_params[c("parallel", "bpparam_class")]
```

Stage 5 records the task definition:

```r
step5$params$parallel_task
# "shared_model_by_metacell"
```

CPU utilization alone is not proof that every stage is parallel. Stage 2, Stage 3 setup, and Stage 5 model-cache construction contain intentional serial sections.

## 8. Serial troubleshooting

A failing stage should first be rerun serially to obtain a deterministic traceback.

Stage 1:

```r
step1_serial <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "RegCompass_steps/01_grn_serial",
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
  parallel = FALSE,
  BPPARAM = FALSE
)
```

Stage 3:

```r
step3_serial <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_meta_modules_serial",
  layer1_args = list(
    local_fastcore = TRUE,
    local_fastcore_args = list(
      solver = "highs",
      parallel = FALSE,
      backend = "serial"
    )
  )
)
```

Stage 5:

```r
step5_serial <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2_serial",
  layer2_args = list(solver = "highs"),
  parallel = FALSE,
  BPPARAM = FALSE
)
```

Never use `BPPARAM = TRUE`; use `FALSE`, `NULL`, or a valid `BiocParallelParam` object.

## 9. Rerun only the medium and Layer 2

```r
medium_scenarios_alt <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "low_glucose",
  species = "human"
)

step5_low_glucose <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios_alt,
  outdir = "RegCompass_steps/05_layer2_low_glucose",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)

result_low_glucose <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5_low_glucose,
  gem = gem,
  outdir = "RegCompass_steps/06_results_low_glucose"
)
```

Do not compare two medium scenarios unless the GEM, target reactions, Layer 1 reaction expression, target directions, solver settings, and convergence criteria are otherwise identical. Worker count may differ because it should not change the objective definition.

## 10. Compare meta-module and full-GEM scoring

`meta_module_gem` scores the shared global union of GRN-derived meta-modules after feasibility completion. `full_gem` keeps the complete prepared GEM and uses the same target core reactions.

```r
step5_full <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2_full_gem",
  model_mode = "full_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 60
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

Interpret differences as sensitivity to structural model scope. Do not treat agreement between the two modes as independent biological replication.

## 11. Solver selection

The default solver is HiGHS and the `highs` R package is a required dependency.

```r
stopifnot(requireNamespace("highs", quietly = TRUE))
```

Optional backends require their matching R packages:

| `solver` value | Required package |
|---|---|
| `"highs"` | `highs` |
| `"glpk"` | `Rglpk` |
| `"gurobi"` | `gurobi` plus a valid Gurobi installation and license |

A solver-installation error is not evidence that the GEM or medium is infeasible. Solver libraries can also start internal threads; keep the BLAS/OpenMP environment controls above when using several R workers.

## 12. Pando installation diagnostics

A local source archive is supported:

```r
install.packages(
  "~/Pando_regcompass.tar.gz",
  repos = NULL,
  type = "source"
)
```

Confirm the required API:

```r
stopifnot(
  requireNamespace("Pando", quietly = TRUE),
  all(vapply(
    c("initiate_grn", "find_motifs", "infer_grn", "gof"),
    exists,
    logical(1),
    envir = asNamespace("Pando"),
    inherits = FALSE
  ))
)
```

GitHub remote metadata are optional. A local installation is accepted when the required API is present.

## 13. Failure classification

| Message or symptom | Classification | Action |
|---|---|---|
| normalized RNA or ATAC contains different cells | input alignment | compare cell-ID sets across assays and metadata |
| `Some features contain 0 total counts` from an old installation | stale package code | reinstall the current source; cell-type-local zero peaks should be excluded before TF-IDF |
| a Pando group is skipped | insufficient cells | inspect `pando_group_status.tsv.gz`; change `min_cells` only with a documented rationale |
| no significant Pando edges | GRN evidence failure | inspect all coefficients, FDR, R², target-gene overlap, and motif inputs |
| dominant cell-type tie | metacell assignment ambiguity | inspect membership/composition; do not assign a GRN arbitrarily |
| GRN/metacell groups do not align | stage contract failure | inspect `grn_metacell_group_coverage.tsv.gz`; rerun the changed upstream stage and Stage 3 |
| no complete-GPR core reaction | mapping failure | inspect GEM gene symbols, GPR parsing, metabolic gene overlap, and projected nodes |
| a worker reports a remote error | parallel execution failure | rerun the same stage serially and inspect the original condition or module |
| worker process is killed without R error | memory or scheduler failure | reduce worker count; inspect kernel, cgroup, or scheduler logs |
| solver package is missing | installation failure | install the selected solver backend |
| parent GEM is `infeasible` after solver preflight | structural or medium failure | inspect applied exchange bounds and parent-model diagnostics |
| some target directions are blocked | possible biological/structural result | inspect target bounds and `lp_diagnostics`; do not convert blocked directions into zeros silently |

## 14. Distinguish medium infeasibility from target blockage

```r
step5$model_cache_summary
head(step5$model_diagnostics)
head(step5$lp_diagnostics)
```

- Parent-model infeasibility means the medium-constrained structural model itself failed.
- Target blockage means the parent model exists but a particular reaction direction cannot carry the required flux.
- A missing solver or solver error means feasibility was not established and must not be interpreted biologically.
- A killed parallel worker is a resource failure until the same task is reproduced serially.

## 15. Output interpretation boundaries

- Pando coefficients are learned from single cells within each condition × cell-type group.
- Metacells are descriptive pseudo-observations built only within condition and guided before aggregation by the selected annotation label.
- The audited dominant member-cell type selects the matching condition × cell-type GRN for Layer 1; inspect composition and purity even when label guidance is enabled.
- Condition contrasts are not biological-sample-level significance tests.
- Sample metadata are provenance only; no sample balancing, weighting, or downsampling is performed.
- FASTCORE support reactions ensure local feasibility and are not additional GRN-supported core reactions.
- Parallelism changes execution order and resource usage, not the biological evidence model.

## 16. Archive a reproducible run

Keep these together:

```text
package version and installation source
Linux distribution and kernel
R version and BiocParallel version
worker counts and backend
OMP/BLAS/MKL environment variables
Seurat object checksum or immutable input path
GEM model_info and checksum
medium_scenarios table
all stage wrapper RDS files
Pando status and edge tables
metacell membership and composition tables
core-reaction and module-membership tables
Layer 2 model cache and diagnostics
final result RDS
```

The one-shot workflow writes `00_model_info.rds` and `00_medium_scenarios.rds` at the run root. Preserve them with the result.
