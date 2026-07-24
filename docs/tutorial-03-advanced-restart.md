# Tutorial Level 3: restart, sensitivity, and diagnostics

Use saved classed stage objects. RegCompassR 1.8.2 rejects cross-run object mixing when GEM fingerprints, workflow parameters, stage classes, core sets, or ordered scoring units differ.

## Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
result <- readRDS("RegCompass_steps/06_results/regcompass_result.rds")
```

Use the stage wrapper RDS, not a compact inspection artifact. Keep all files referenced by `step5$model_cache_summary$file`.

## Earliest stage to rerun

| Change | Rerun from |
|---|---:|
| Pando model, `tf_cor`, `peak_cor`, minimum cells | Stage 1 |
| metacell `gamma` or cell-size thresholds | Stage 2 |
| condition/cell-type metadata or assay names | Stage 1 and Stage 2 |
| GRN projection, subsystem/cross-reference expansion, local FASTCORE | Stage 3 |
| regulatory integration, GPR `tau`, RNA half-saturation | Stage 4 |
| medium, structural mode, solver, `omega`, target direction | Stage 5 |
| selected core anchors for direct database-linked scoring | target-union step only |
| ranking/annotation assembly | Stage 6 |

A changed GEM invalidates Stage 1, Stage 3, Stage 4, Stage 5, Stage 6, and target-union outputs because their fingerprints no longer match.

## Linux workers

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export REGCOMPASS_WORKERS=16
```

```r
library(BiocParallel)
upstream_bp <- MulticoreParam(workers = 16L, progressbar = TRUE)
layer2_bp <- MulticoreParam(workers = 12L, progressbar = TRUE)
```

| Stage | Parallel unit |
|---|---|
| 1 | condition × cell-type Pando group |
| 2 | no workflow-level BiocParallel loop |
| 3 | local FASTCORE completion per meta-module |
| 4 | GPR/reaction-capacity calculation |
| 5 | shared model × metacell |
| target union | reused union model × metacell |
| 6 | serial assembly |

Keep Pando's inner `parallel = FALSE`. Lower Layer 2 worker counts when memory, not CPU, is limiting.

## Restart Stage 5 with a new medium

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = c("normal_human_plasma", "dmem_high_glucose"),
  species = "human"
)

step5_medium <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/05_layer2_medium",
  model_mode = "meta_module_gem",
  layer2_args = list(solver = "highs", target_direction = "both"),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

Medium bounds can restrict existing GEM directions but cannot create a direction absent from the source model.

## Restart direct database-linked scoring

```r
expanded <- rc_regcompass_step_target_union(
  layer1 = step4,
  meta_modules = step3,
  layer2 = step5,
  gem = gem,
  outdir = "RegCompass_steps/05b_glutathione",
  core_reaction_ids = "MAR04324",
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

The selected core is not scored again. Only non-core reactions directly sharing a KEGG, Reactome, or master-Rhea ID with the selected core are evaluated. Same-subsystem and transitive expansion are unavailable by design.

This restart is valid only while the original union-GEM cache files remain unchanged and available. The output records their paths, fingerprints, MD5 hashes, and sizes.

## Serial troubleshooting

Rerun only the failing computational stage with `parallel = FALSE` and `BPPARAM = FALSE`. For Stage 3, set `local_fastcore_args$parallel = FALSE` and `backend = "serial"`.

Classify failures in this order:

1. **Input contract:** missing assays, metadata, stage class, fingerprint, or reordered units.
2. **Installation:** unavailable Pando, SuperCell2, genome package, or LP solver.
3. **Model construction:** missing core reactions, invalid GPRs, or incomplete FASTCORE support.
4. **Database mapping:** selected cores have no direct KEGG, Reactome, or master-Rhea-linked non-core reactions in the original union.
5. **Medium:** exchange mapping or restrictive bounds.
6. **Target direction:** reaction blocked in the fixed model.
7. **Evidence:** finite LP result but weak or non-informative variation.

Do not interpret a solver-installation error as biological infeasibility. Do not interpret a blocked direction as evidence that the opposite direction is inactive without checking the signed model bounds.

## Validate a restart before interpretation

```r
stopifnot(
  inherits(step3, "regcompass_meta_module_step"),
  inherits(step4, "regcompass_layer1_step"),
  inherits(step5, "regcompass_layer2_step"),
  identical(step3$gem_fingerprint, step4$gem_fingerprint),
  identical(step4$gem_fingerprint, step5$gem_fingerprint),
  identical(step4$workflow_params, step5$workflow_params),
  identical(colnames(step4$reaction_expression), colnames(step5$penalty)),
  all(file.exists(step5$model_cache_summary$file))
)
```

See [Predefined extracellular medium scenarios](medium-presets.md) and [Direct database-linked non-core scoring](target-union-scoring.md) for the corresponding contracts.
