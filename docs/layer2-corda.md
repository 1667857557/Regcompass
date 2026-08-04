# Optional original CORDA reconstruction in Layer 2

`rc_regcompass_step_layer2()` keeps compact add-only FASTCORE as the default.
For `model_mode = "meta_module_gem"`, `model_completion = "corda"` runs the
three-stage CORDA algorithm described by Schultz and Qutub (2016). The former
draft name `"corda_like"` is accepted as an alias and normalized to `"corda"`.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_corda",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "corda",
      completion_time_limit = 3000,
      fastcore_epsilon = 1e-4,
      strict = TRUE,
      corda_gamma = 1e5,
      corda_kappa = 1e-2,
      corda_epsilon = 1,
      corda_n = 5L,
      corda_p = 2L,
      corda_seed = 1L,
      corda_flux_tolerance = 1e-8,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

## Confidence mapping

CORDA accepts reaction confidence classes; it does not prescribe how omics data
must be converted into those classes. RegCompass performs that mapping within
each cell type.

RNA-only and multiome reaction capacities are summarized across matching
metacells and converted to within-cell-type percentile ranks. The evidence score
is

\[
E_r=(1-w)\max\{Q_r^{RNA},Q_r^{multiome}\}+wA_r^{regulatory},
\]

where `A_regulatory` is the median GPR regulatory-support fraction and the
default `w` is `0.20`.

The initial classes are:

- **HC:** merged core reactions;
- **MC_module:** non-core reactions in the merged cell-type meta-module;
- **MC_evidence:** optional reactions outside the module with evidence at or
  above `corda_medium_confidence_threshold`;
- **NC:** remaining reactions with finite evidence at or below
  `corda_negative_confidence_threshold`;
- **OT:** all remaining reactions, including reactions with unavailable omics
  evidence.

MC reactions are flexible. They are not automatically retained.

## Directional representation

The medium-constrained, FASTCC-audited parent GEM is split once into nonnegative
forward and reverse variables. The original reaction bounds determine which
directions exist and their directional upper and lower bounds. CORDA operates on
these directional variables, while the final model is converted back to original
reaction IDs, original stoichiometry, bounds, annotations and GPR rules.

## Dependency assessment

For a target direction, CORDA requires at least `corda_epsilon` flux and solves

\[
\min_v \sum_j c_jv_j,
\qquad S_{split}v=0,
\qquad l\le v\le u.
\]

Undesirable directions receive their stage-specific base cost. Every direction
also receives deterministic uniform noise in `[0, corda_kappa]`. The noise is
keyed by stage, target, direction, repeat and `corda_seed`, so results do not
depend on worker scheduling. Penalized reactions carrying flux above
`corda_flux_tolerance` are associated with the target.

The paper defaults are:

- `corda_gamma = 1e5`;
- `corda_kappa = 1e-2`;
- `corda_epsilon = 1`;
- `corda_n = 5` randomized dependency assessments per target direction;
- `corda_p = 2` distinct MC reactions required to promote a shared NC reaction.

## Three reconstruction stages

### Stage 1: HC dependencies

All HC directions are assessed `n` times. MC directions have cost
`sqrt(gamma)`, NC directions have cost `gamma`, and HC/OT directions have zero
base cost. Any MC or NC reaction associated with an HC target is promoted to the
retained set.

### Stage 2: flexible MC and shared NC support

Each remaining MC direction is assessed `n` times while only NC directions are
penalized by `gamma`. An NC reaction is promoted when it is associated with at
least `p` distinct MC reactions. Remaining NC reactions are then blocked. Every
remaining MC direction is tested for flux capacity, and an MC reaction is
promoted when at least one allowed direction can carry `corda_epsilon` flux.

### Stage 3: OT dependencies of the retained set

All still-unpromoted MC and NC reactions are blocked. Every retained direction
is assessed `n` times while OT directions have cost `gamma`. Associated OT
reactions are promoted. The final model contains the retained reaction set.

Confidence updates occur only after all tasks in a stage finish. This barrier
semantics makes the reconstruction independent of task completion order.

## Native C++ acceleration and parallelism

The LP solve is the dominant cost. The `highs` package already exposes the
native HiGHS C++ solver. RegCompass therefore does not duplicate the solver in a
custom Rcpp implementation. Instead, each worker:

1. constructs one persistent HiGHS model for its task block;
2. updates only the objective coefficients and variable bounds between tasks;
3. re-solves with the simplex solver and reuses the existing basis;
4. falls back to the previous one-shot LP path if the installed `highs` version
   does not expose the persistent solver API or a persistent call fails.

When the number of cell-type-by-medium models is at least the worker count,
models run in parallel and each model executes its stages serially inside one
worker. When there are fewer models than workers, models run sequentially and
each stage distributes `target direction × repeat` chunks across the full worker
pool. A supplied but inactive `BPPARAM` is started once for the whole Layer 2
call and released at exit; an already active caller pool remains active. When
`BPPARAM = NULL`, the CORDA route creates one package-managed pool and keeps it
alive across all stages.

The persistent native solver removes repeated model construction and permits
basis reuse, but the actual speedup depends on model size, presolve behavior and
worker count. Performance claims should be based on a Human-GEM benchmark rather
than the synthetic correctness test.

## Pipeline contracts

- Only parent-feasible HC core directions enter the downstream `vmax`, penalty
  and score matrices. MC/NC/OT reconstruction does not change matrix row names.
- Conditions and matching metacells share one model only within the same cell
  type and medium.
- The RNA-only control reuses the exact CORDA model cache and checksum.
- Final models retain original reaction IDs, bounds, annotations and GPR rules.
- CORDA caches are isolated under
  `model_cache/meta_module_gem/corda`.
- `max_support_reactions` remains a FASTCORE control and is recorded but not used
  to truncate CORDA, because such a cap would change the published algorithm.

## Audit output

Each model records initial confidence, final retained status, inclusion stage,
reaction evidence, dependency tasks, Stage 2 NC-to-MC support pairs, execution
backend, persistent-solver fallback counts and core-direction closure checks.
The full task table is stored once and a compact stage/status/backend summary is
stored separately to avoid duplicate model-cache payloads. The Layer 2 result
records `params$model_completion` and `completion_contract`; the exported
model-cache summary reports class and stage counts for every cell-type-by-medium
model.
