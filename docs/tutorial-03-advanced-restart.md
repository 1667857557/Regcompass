# Tutorial 3: restart and sensitivity

RegCompass saves one object per stage. Restart from the earliest stage affected
by a changed input or parameter.

## Load a completed run

```r
step1 <- readRDS("RegCompass_steps/01_grn/step_grn.rds")
step2 <- readRDS("RegCompass_steps/02_metacells/step_metacells.rds")
step3 <- readRDS("RegCompass_steps/03_meta_modules/step_meta_modules.rds")
step4 <- readRDS("RegCompass_steps/04_layer1/step_layer1.rds")
step5 <- readRDS("RegCompass_steps/05_layer2/step_layer2.rds")
```

## Restart boundaries

| Restart from | Changes |
|---|---|
| Stage 1 | RNA/ATAC data, condition or cell-type labels, genome, regulatory regions, motifs, Pando fitting or filtering arguments |
| Stage 2 | reductions, dimensions, gamma, seed, metacell thresholds, or cells |
| Stage 3 | GPR rules, subsystem table, KEGG/Reactome mappings, or master-Rhea mappings |
| Stage 4 | `regulatory_alpha`, `gpr_and_method`, gene half-saturation, or metacell evidence |
| Stage 5 | medium, exchange bounds, target direction, omega, solver, or model-completion settings |
| Stage 6 | annotations, reporting filters, or contrast settings only |

Changing an earlier stage invalidates every downstream stage.

## Stage 4 sensitivity

The default GPR AND rule is `"min"`. Use `"median"` or `"mean"` to test
sensitivity to complex aggregation.

```r
step4_mean <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "RegCompass_restart/04_layer1_mean",
  projection_component = "condition",
  comparison_support = "auto",
  regulatory_alpha = 0.5,
  gpr_and_method = "mean",
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Other useful Stage 4 sensitivity settings:

- `regulatory_alpha`: regulatory contribution to RNA support;
- `gene_half_saturation`: RNA support saturation;
- `comparison_support`: `pairwise_common` or `global_common` when explicit control is required.

## Change medium

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
  outdir = "RegCompass_restart/05_layer2_low_glucose",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
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

Medium choices:

- `physiologic`: species-specific baseline;
- `normal_human_plasma` or `mouse_plasma`: explicit physiological preset;
- `rpmi1640` or `dmem_high_glucose`: culture formulation;
- nutrient challenge presets: human-only sensitivity scenarios;
- `minimal`, `compass_model_bounds`, `permissive_all_exchange`: technical sensitivity scenarios;
- `custom`: experiment-specific medium.

See [Medium presets](medium-presets.md) for species restrictions and preset
contents.

## Model-completion diagnostics

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
table(model$closure_diagnostics$completion_status)
```

Common statuses:

- `already_feasible`: no additional support reaction was required;
- `global_fastcore_completed`: FASTCORE added support reactions;
- `parent_blocked`: target direction is infeasible in the parent model;
- `unresolved`: completion failed under the configured limits;
- `no_allowed_direction`: original bounds block the direction.

`completion_time_limit` applies to model completion, not the later scoring LPs.

## Reassemble results

```r
result_new <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5_new,
  gem = gem,
  outdir = "RegCompass_restart/06_results_low_glucose",
  species = "human"
)
```

Mathematical definitions: [Mathematical model](mathematical-model.md). Public
API: [functions.md](functions.md).
