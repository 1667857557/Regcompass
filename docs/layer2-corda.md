# Original CORDA2 in Layer 2

`rc_regcompass_step_layer2()` uses original MATLAB CORDA2 semantics by default when `model_mode = "meta_module_gem"`. Set `model_params$model_completion = "fastcore"` only when the supplementary FASTCORE route is required.

## Parameters

```r
layer2_args = list(
  target_direction = "both",
  solver = "highs",
  model_params = list(
    completion_time_limit = 3000,
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
```

| Argument | Default | Meaning |
|---|---:|---|
| `MCxNCthresh` | `2` | Minimum number of MC directional dependencies required to promote an NC direction. |
| `constraint` | `1` | Absolute target flux for `constrainby = "val"`, or percentage for `"perc"`. |
| `constrainby` | `"val"` | Target constraint mode. |
| `om` | `1e4` | High reaction cost; Step-1 MC cost is `sqrt(om)`. |
| `ci` | `0.01` | Proportional cost increase for newly observed high-cost dependencies. |

## RegCompass confidence adapter

For each `reaction × cell type`, Layer 1 evidence is summarized across metacells and conditions and mapped to the four CORDA2 groups:

| RegCompass evidence class | CORDA2 group |
|---|---|
| merged core reaction | HC |
| merged non-core or selected high-evidence reaction | MC |
| finite low-evidence reaction | NC |
| remaining reaction | OT |

The reconstruction then follows the original three-stage CORDA2 state machine. The complete medium-constrained parent GEM is passed directly to CORDA2 without FASTCC or role-based pre-pruning.

## Directional bounds

Only reactions with `lb < 0 && ub >= 0` are split into forward and `_CORDA_rev_rxn` variables. The opposite direction is closed while a target direction is assessed.

After reconstruction, selected directional variables are merged back to the original reaction identifiers. The final model restores the medium-constrained parent bounds for retained directions:

- a retained reversible direction keeps the corresponding parent directional bound;
- a retained irreversible reaction with `parent lb > 0` recovers that positive lower bound;
- an excluded direction is fixed to zero.

This avoids converting a required positive parent flux into an optional `lb = 0` reaction during finalization.

## Penalty scale

Layer 2 scoring uses the COMPASS cost function:

```text
P = 1 / (1 + log2(1 + max(E, 0)))
```

Missing expression and structural roles use the maximum COMPASS cost `1`; no structural role receives a default cost of `20`.

## Execution and provenance

One CORDA2 reconstruction is built per `cell type × medium`. Independent models may run through `BPPARAM`; target assessments within one reconstruction remain serial. The completed model is reused for primary multiome and RNA-only scoring.

Inspect `step5$completion_contract`, `step5$model_cache_summary`, `step5$vmax_cache_diagnostics` and `step5$lp_diagnostics` for the algorithm, parameters and target feasibility results.
