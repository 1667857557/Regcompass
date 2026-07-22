# Tutorial Level 3: restart, sensitivity runs, and diagnostics

Use this level after completing the [Level 2 audit workflow](tutorial-02-stepwise-audit.md). It explains which stages must be rerun after a parameter change and how to separate installation, structural, and biological failures.

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

Always load the stage wrapper, such as `step_grn.rds` or `step_metacells.rds`, when a downstream function requires the class and workflow parameters. The compact files such as `single_cell_grn.rds` are useful for inspection but are not substitutes for the stage wrapper.

## 2. Minimal rerun matrix

| Changed item | Earliest stage to rerun | Reusable upstream stages |
|---|---:|---|
| Pando `tf_cor`, `peak_cor`, model method, minimum cells | 1 | none for GRN-derived modules; Stage 2 may be reused if unchanged |
| metacell `gamma`, minimum stratum size, minimum metacell size | 2 | Stage 1 |
| GRN projection or meta-module expansion settings | 3 | Stages 1-2 |
| `regulatory_alpha`, GPR `tau`, RNA half-saturation | 4 | Stages 1-3 |
| medium scenario, solver, `omega`, target direction, `model_mode` | 5 | Stages 1-4 |
| display or final assembly only | 6 | Stages 1-5 |

When Stage 1 or Stage 2 changes, rerun Stage 3 because the bidirectional GRN/metacell group-coverage contract must be revalidated.

## 3. Rerun only the medium and Layer 2

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
  )
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

Do not compare two medium scenarios unless the GEM, target reactions, Layer 1 reaction expression, target directions, and solver settings are otherwise identical.

## 4. Compare meta-module and full-GEM scoring

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
  )
)
```

Interpret differences as sensitivity to structural model scope. Do not treat agreement between the two modes as independent biological replication.

## 5. Parallel execution

For reproducible troubleshooting, first rerun a failing stage serially:

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
      adjust_method = "fdr"
    )
  ),
  parallel = FALSE,
  BPPARAM = FALSE
)
```

For the one-shot workflow, resource controls are exposed through `upstream_workers`, `layer2_workers`, and `parallel_backend`. Do not use a logical value such as `BPPARAM = TRUE`; use `FALSE`, `NULL`, or a valid BiocParallel parameter object.

## 6. Solver selection

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

A solver-installation error is not evidence that the GEM or medium is infeasible.

## 7. Pando installation diagnostics

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

## 8. Failure classification

| Message or symptom | Classification | Action |
|---|---|---|
| normalized RNA or ATAC contains different cells | input alignment | compare cell-ID sets across assays and metadata |
| `Some features contain 0 total counts` from an old installation | stale package code | reinstall the current source; cell-type-local zero peaks should be excluded before TF-IDF |
| a Pando group is skipped | insufficient cells | inspect `pando_group_status.tsv.gz`; change `min_cells` only with a documented rationale |
| no significant Pando edges | GRN evidence failure | inspect all coefficients, FDR, R², target-gene overlap, and motif inputs |
| dominant cell-type tie | metacell assignment ambiguity | inspect membership/composition; do not assign a GRN arbitrarily |
| GRN/metacell groups do not align | stage contract failure | inspect `grn_metacell_group_coverage.tsv.gz`; rerun the changed upstream stage and Stage 3 |
| no complete-GPR core reaction | mapping failure | inspect GEM gene symbols, GPR parsing, metabolic gene overlap, and projected nodes |
| solver package is missing | installation failure | install the selected solver backend |
| parent GEM is `infeasible` after solver preflight | structural or medium failure | inspect applied exchange bounds and parent-model diagnostics |
| some target directions are blocked | possible biological/structural result | inspect target bounds and `lp_diagnostics`; do not convert blocked directions into zeros silently |

## 9. Distinguish medium infeasibility from target blockage

```r
step5$model_cache_summary
head(step5$model_diagnostics)
head(step5$lp_diagnostics)
```

- Parent-model infeasibility means the medium-constrained structural model itself failed.
- Target blockage means the parent model exists but a particular reaction direction cannot carry the required flux.
- A missing solver or solver error means feasibility was not established and must not be interpreted biologically.

## 10. Output interpretation boundaries

- Pando coefficients are learned from single cells within each condition × cell-type group.
- Metacells are descriptive pseudo-observations built only within condition.
- The post hoc dominant cell type selects the matching condition × cell-type GRN for Layer 1.
- Condition contrasts are not biological-sample-level significance tests.
- Sample metadata are provenance only; no sample balancing, weighting, or downsampling is performed.
- FASTCORE support reactions ensure local feasibility and are not additional GRN-supported core reactions.

## 11. Archive a reproducible run

Keep these together:

```text
package version and installation source
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
