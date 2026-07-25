# Tutorial Level 3: restart, sensitivity, and diagnostics

RegCompass saves classed stage objects so downstream stages can be rerun without repeating unchanged work. Objects from incompatible runs are rejected through stored workflow and GEM provenance.

## Load a completed stepwise run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

## Restart boundaries

### Rerun Stage 2 onward

Rerun metacells and all downstream stages after changing cell membership, RNA/ATAC matrices, condition or cell-type metadata, `gamma`, or the RNA/ATAC reductions used for metacell construction.

### Rerun Stage 3 onward

Rerun reaction meta-modules and downstream stages after changing Pando filtering, GRN-neighbour parameters, expansion settings, subsystem/cross-reference annotations, or GEM GPR rules.

### Rerun Stage 4 onward

Rerun Layer 1 and downstream stages after changing `regulatory_alpha`, `tau`, gene half-saturation, or metacell RNA/ATAC evidence while retaining the same Stage 3 reaction catalogue.

### Rerun Stage 5 onward

Rerun Layer 2 and results after changing:

- medium composition or exchange bounds;
- `target_direction`, `omega`, solver, or scoring `time_limit`;
- `completion_time_limit`, `fastcore_epsilon`, `max_support_reactions`, or `strict` in `layer2_args$model_params`.

The complete preset list and custom-medium format are documented in [medium presets](medium-presets.md).

## Rebuild Stage 5

```r
new_medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "low_glucose",
  species = "human"
)

step5_new <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = new_medium_scenarios,
  outdir = "RegCompass_restart/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    time_limit = 900,
    model_params = list(
      completion_time_limit = 900,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2500,
      strict = TRUE
    )
  ),
  parallel = TRUE,
  BPPARAM = layer2_bp
)
```

This reuses the Stage 3 reaction targets and Stage 4 support matrix.

## Diagnose model completion

```r
summary <- step5_new$model_cache_summary
summary[, c(
  "medium_scenario",
  "n_merged_biological_reactions",
  "n_global_fastcore_support_reactions",
  "target_status",
  "file",
  "file_checksum"
)]

model <- readRDS(summary$file[[1]])
diagnostics <- model$closure_diagnostics
table(diagnostics$completion_status)
```

Common statuses are:

- `already_feasible`: the initial reaction set supports the target;
- `global_fastcore_completed`: FASTCORE added supporting reactions;
- `parent_blocked`: the target direction is infeasible in the medium-constrained parent GEM;
- `unresolved`: completion did not succeed under the requested limits;
- `no_allowed_direction`: the original GEM bounds block the requested direction.
