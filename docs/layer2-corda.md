# Original CORDA2 in Layer 2

`rc_regcompass_step_layer2()` uses original MATLAB CORDA2 semantics by default when `model_mode = "meta_module_gem"`. Set `model_params$model_completion = "fastcore"` only when the supplementary FASTCORE route is required.

## Parallel worker cap

Layer 2 uses the same single RegCompass worker parameter as the rest of the workflow:

```r
workers = 10L
```

The default requested cap is `10L` and users may change it. RegCompass automatically chooses `BiocParallel::SnowParam(type = "SOCK")` on Windows and `BiocParallel::MulticoreParam` on Linux/macOS. Two detected logical CPUs are reserved, so the effective cap is

```text
min(workers, max(1, detected logical CPUs - 2))
```

Each CORDA2 step additionally uses no more workers than it has independent directional targets.

## Parameters

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  workers = 10L,
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

CORDA2 reconstruction intentionally has no structural time limit. `model_params$completion_time_limit` is rejected on the CORDA2 route so a long reconstruction cannot be silently truncated. That parameter remains available only for supplementary non-CORDA2 completion such as FASTCORE.

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

The reconstruction follows the original CORDA2 state machine. The complete medium-constrained parent GEM is passed directly to CORDA2 without FASTCC or role-based pre-pruning.

## Stage-barrier target parallelism

The mathematical steps remain strictly sequential:

```text
Step 1 HC dependencies
  -> ordered reduction and HC/MC/NC update
Step 2.1 MC-to-NC dependencies
  -> ordered reduction and MCxNC update
Step 2.2 MC feasibility
  -> ordered reduction and NC/MC promotion update
Step 3 HC-to-OT dependencies
  -> ordered reduction and final OT inclusion
```

Within one step, directional targets are independent with respect to that step's frozen pre-step state and are distributed across the protected worker budget. Every worker uses the unchanged CORDA2 dependency/maximization routines and one HiGHS thread. Results are restored to the original directional-target indices before any dependency matrix, promotion, blocking or confidence state is modified.

After every step RegCompass releases worker-local HiGHS engines, stops the step-local worker pool, performs full garbage collection and only then enters the next step. Independent cell-type × medium CORDA2 models are processed one at a time so the current model can use the full protected worker budget inside its mathematical step.

## Progress

Entering each step prints its name, target count, worker count and chunk count. During the step the progress display reports

```text
completed / total
percentage
remaining directional targets
current completed target
```

Step completion reports `remaining=0` and confirms that the worker pool was released. The same events are persisted through the Layer 2 progress infrastructure with `scope = "corda2_stage"`.

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

One CORDA2 reconstruction is built per `cell type × medium` and reused for primary multiome and RNA-only scoring. Inspect `step5$completion_contract`, `step5$model_cache_summary`, `step5$vmax_cache_diagnostics`, `step5$lp_diagnostics` and each reconstructed model's `corda_reconstruction$stage_parallelism` for the algorithm, worker counts, stage target counts, parameters and target-feasibility results.
