# Layer 2 structural model builders

`rc_regcompass_step_layer2()` supports three mutually exclusive structural routes. In every route, a medium scenario changes listed exchange-reaction bounds only; the same completed structural cache is reused for the primary multiome score and the RNA-only control.

## Default: CORDA2 cell-type models

With `model_mode = "meta_module_gem"`, omitting `model_params$model_completion` selects `"corda2"`.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      strict = TRUE,
      corda2_args = list(
        MCxNCthresh = 2,
        constraint = 1,
        constrainby = "val",
        om = 1e4,
        ci = 0.01
      ),
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

The complete medium-constrained parent is passed directly to the original MATLAB CORDA2 state machine. FASTCC and role-based pre-pruning are not run. One model is reconstructed per `cell type × medium` and shared across conditions of that cell type. During directional merge, selected reactions recover their parent bounds, including positive lower bounds.

CORDA2 reconstruction is intentionally unbounded in wall-clock time. `model_params$completion_time_limit` is rejected when CORDA2 is selected so a large reconstruction cannot be partially completed and mistaken for a valid context-specific model.

## Supplementary: FASTCORE

```r
step5_fastcore <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_fastcore",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "fastcore",
      completion_time_limit = 1200,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 3000,
      strict = TRUE
    )
  )
)
```

FASTCORE is now an explicit alternative. It performs FASTCC consistency analysis and add-only support selection after medium bounds are applied. `completion_time_limit` remains available for this supplementary non-CORDA2 route.

## Supplementary: COMPASS-style full GEM

```r
step5_full <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_full_gem",
  model_mode = "full_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    flux_threshold = 1e-8
  )
)
```

Full-GEM mode keeps the complete reference GEM and skips FASTCC, FASTCORE and CORDA2. Do not supply completion-specific controls in this mode.

## Penalty and output contract

All routes use the COMPASS reaction-cost scale:

```text
P = 1 / (1 + log2(1 + max(E, 0)))
```

Missing reaction expression and the structural roles `exchange`, `demand`, `sink` and `artificial_support` receive the maximum cost `1`.

| Route | `model_completion` | default | `fastcore_executed` | `corda2_executed` | structural time limit |
|---|---|---:|---:|---:|---|
| CORDA2 | `"corda2"` | yes | `FALSE` | `TRUE` | none |
| FASTCORE | `"fastcore"` | no | `TRUE` | `FALSE` | optional `completion_time_limit` |
| Full GEM | `"none"` | no | `FALSE` | `FALSE` | not applicable |

Inspect `step5$params`, `step5$completion_contract`, `step5$model_cache_summary`, `step5$vmax_cache_diagnostics` and `step5$lp_diagnostics` for the selected route and feasibility results.
