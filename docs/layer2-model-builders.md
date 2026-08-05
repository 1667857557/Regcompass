# Layer 2 structural model builders

`rc_regcompass_step_layer2()` supports exactly three structural choices. They are mutually exclusive and all reuse the same completed structural cache for the primary multiome score and the RNA-only control.

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

Omitting `model_completion` under `model_mode = "meta_module_gem"` is equivalent to `model_completion = "fastcore"`. One cell-type-by-medium biological reaction union is completed by add-only FASTCORE. FASTCC is used to define the consistent parent model before support reactions are selected.

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

CORDA2 is available only with `model_mode = "meta_module_gem"`. FASTCORE is not executed on this route. Independent cell-type-by-medium CORDA2 reconstructions may run in parallel; target processing inside each reconstruction remains serial.

## 3. Medium-pruned full GEM

The full-GEM route applies the requested medium to the reference GEM and removes reactions that cannot carry non-zero steady-state flux under that medium. It does not use reaction-expression evidence to construct the model and does not execute FASTCORE or CORDA2.

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
    flux_threshold = 1e-8,
    model_params = list(
      completion_time_limit = 1200
    )
  )
)
```

`model_completion` should normally be omitted. The explicit value `"none"` is accepted, but `"fastcore"` and `"corda2"` are rejected. FASTCORE- and CORDA2-specific controls are also rejected rather than silently ignored.

The structural sequence is:

```text
reference GEM
  -> apply shared medium bounds
  -> FASTCC flux-consistency analysis
  -> remove reactions with no feasible non-zero steady-state flux
  -> remove orphan metabolites
  -> directional COMPASS-like scoring
```

FASTCC here is only a consistency filter. It does not select a compact evidence-supported network and therefore is not a FASTCORE reconstruction. The workflow still uses the Stage 3 merged core reactions as scoring targets; the pruned full GEM supplies the complete feasible support network for those targets.

## Output contract

The selected route is recorded in `step5$params` and `step5$completion_contract`:

| Route | `model_completion` | `fastcore_executed` | `corda2_executed` |
|---|---|---:|---:|
| FASTCORE | `"fastcore"` | `TRUE` | `FALSE` |
| CORDA2 | `"corda2"` | `FALSE` | `TRUE` |
| medium-pruned full GEM | `"none"` | `FALSE` | `FALSE` |

For the full-GEM route, `model_cache_summary` also records the input reaction count, retained reaction count, number of removed flux-inconsistent reactions, solver, structural time limit, and consistency threshold. Cache fingerprints include the reference GEM, pruning algorithm, solver, and threshold so an older unpruned cache cannot be reused silently.
