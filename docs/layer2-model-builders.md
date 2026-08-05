# Layer 2 structural model builders

`rc_regcompass_step_layer2()` supports exactly three mutually exclusive structural choices. All three use the same medium contract:

```text
medium scenario
  -> modify bounds of the listed exchange reactions
  -> retain the reaction set and stoichiometric matrix
  -> let the selected algorithm or directional LP determine feasibility
```

The medium table itself never removes reaction or metabolite columns. The completed structural cache is reused for both the primary multiome score and the RNA-only control.

## 1. FASTCORE

FASTCORE is the default cell-type-specific reconstruction route.

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
      model_completion = "fastcore",
      completion_time_limit = 1200,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 3000,
      strict = TRUE
    )
  )
)
```

Omitting `model_completion` under `model_mode = "meta_module_gem"` is equivalent to `model_completion = "fastcore"`.

The sequence is:

```text
complete reference GEM
  -> apply medium exchange bounds without deleting reactions
  -> FASTCC consistency analysis inside the FASTCORE route
  -> add-only FASTCORE support selection
  -> compact cell-type-by-medium union GEM
```

FASTCC is an internal part of FASTCORE parent preparation. Reactions absent from the final compact model are excluded by FASTCORE reconstruction, not deleted directly because they were absent from the medium table.

## 2. CORDA2

CORDA2 uses the original MATLAB CORDA2 state-machine semantics and the RegCompass confidence adapter.

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
      model_completion = "corda2",
      completion_time_limit = 1200,
      corda2_args = list(
        MCxNCthresh = 2,
        constraint = 1,
        constrainby = "val",
        om = 1e4,
        ci = 0.01
      )
    )
  )
)
```

The sequence is:

```text
complete reference GEM
  -> apply medium exchange bounds without deleting reactions
  -> pass the complete medium-constrained parent directly to CORDA2
  -> original CORDA2 confidence-state reconstruction
```

There is no FASTCC or role-based parent pre-pruning before CORDA2. CORDA2 is available only with `model_mode = "meta_module_gem"`; FASTCORE is not executed on this route. Independent cell-type-by-medium reconstructions may run in parallel, while target processing inside each CORDA2 reconstruction remains serial.

## 3. COMPASS-style full GEM

The full-GEM route retains the complete reference network. It applies the medium only to exchange-reaction bounds and then follows the COMPASS target-scoring pattern.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "full_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    flux_threshold = 1e-8
  )
)
```

`model_completion` should normally be omitted. The explicit value `"none"` is accepted, but `"fastcore"` and `"corda2"` are rejected. FASTCORE- and CORDA2-specific controls are rejected rather than silently ignored.

The sequence is:

```text
complete reference GEM
  -> apply medium exchange bounds without deleting reactions
  -> retain every requested target direction in the cache
  -> maximize each target direction under steady-state constraints
  -> if vmax < flux_threshold: mark direction infeasible and skip Step 2
  -> otherwise constrain target flux to omega * vmax and minimize penalty
```

Full-GEM mode does not run FASTCC, FASTCORE, or CORDA2. A reaction that cannot operate under a medium remains in the stoichiometric matrix; its infeasibility is represented by its directional `vmax` result and LP diagnostics.

## Output contract

The selected route is recorded in `step5$params` and `step5$completion_contract`:

| Route | `model_completion` | medium directly deletes reactions | `fastcore_executed` | `corda2_executed` |
|---|---|---:|---:|---:|
| FASTCORE | `"fastcore"` | `FALSE` | `TRUE` | `FALSE` |
| CORDA2 | `"corda2"` | `FALSE` | `FALSE` | `TRUE` |
| COMPASS-style full GEM | `"none"` | `FALSE` | `FALSE` | `FALSE` |

For full-GEM mode, `model_cache_summary` records the input and retained reaction counts, the invariant zero count of medium-removed reactions, the number of exchange-bound changes, and the exact medium fingerprint. Cache keys include the reference GEM and exact exchange bounds, so a same-named custom medium with changed bounds cannot reuse a stale model.

The target-level diagnostics are stored in `vmax_cache_diagnostics` and `lp_diagnostics`. Medium-infeasible directions remain visible with `feasible = FALSE` and `step2_status = "not_run"` rather than disappearing from the output.
