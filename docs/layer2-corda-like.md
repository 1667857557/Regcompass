# Optional CORDA-like Layer 2 model completion

`rc_regcompass_step_layer2()` keeps the existing compact FASTCORE completion as
the default. An optional evidence-maximizing CORDA-like route is available only
for `model_mode = "meta_module_gem"`.

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
      model_completion = "corda_like",
      completion_time_limit = 3000,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 3000,
      strict = TRUE,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_other_penalty = 1,
      corda_negative_penalty = 10,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

## Evidence classes

For each cell type, reaction evidence is calculated from the matching Layer 1
metacells. RNA-only and multiome reaction capacities are summarized by their
median and converted to within-cell-type percentile ranks. The default evidence
score is

\[
E_r = 0.8\max\{Q_r^{RNA},Q_r^{multiome}\}
      +0.2A_r^{regulatory},
\]

where `A_regulatory` is the median GPR regulatory-support fraction across the
matching metacells. The `corda_regulatory_weight` parameter replaces `0.2` and
the expression component receives weight `1 - corda_regulatory_weight`.

The structural classes are:

- **HC:** merged core reactions. All HC reactions are retained and their allowed
  directions are checked explicitly.
- **MC_module:** all non-core reactions in the cell-type merged meta-module.
- **MC_evidence:** reactions outside the merged meta-module whose evidence score
  reaches `corda_medium_confidence_threshold`. These reactions are optional as a
  class, but every selected MC reaction is retained in the resulting model.
- **OT:** remaining parent-GEM reactions with missing or intermediate evidence.
- **NC:** remaining reactions with finite evidence at or below
  `corda_negative_confidence_threshold`.

The medium-constrained parent GEM is FASTCC-audited and inconsistent reactions
are fixed to zero bounds. The initial cell-type model retains all HC, MC_module
and MC_evidence entries. Weighted directional completion then restores every
parent-feasible HC/MC direction; only HC directions remain downstream scoring
targets, so expanding MC structure does not enlarge the Layer 2 score matrix.

## Weighted support completion

This implementation is CORDA-like rather than a verbatim implementation of the
original CORDA MILP. It uses the existing add-only LP7/LP10 completion framework,
but replaces the uniform support objective with a weighted L1 objective:

\[
\min \sum_{r\in OT} c_{OT}|v_r|
     +\sum_{r\in NC} c_{NC}|v_r|,
\qquad c_{NC}\ge c_{OT}\ge0.
\]

By default, `c_OT = 1` and `c_NC = 10`. Thus all evidence-supported HC/MC
reactions are retained and made directionally flux-consistent whenever that
direction is feasible in the audited parent model. NC support is strongly
disfavored relative to OT support. Because the objective is weighted flux rather
than a binary reaction count, an NC reaction can still be selected when required
for feasibility or when it gives the lower total weighted support flux.

This objective does not claim that every retained reaction is active in one
common flux vector. Instead, each retained HC/MC direction is required to carry
at least `fastcore_epsilon` in at least one admissible steady state.

## Backward compatibility

Omitting `model_completion`, or setting it to `"fastcore"`, calls the previous
FASTCORE cache builder unchanged. The public function signature, output matrix
ordering, shared-model reuse, directional `vmax` cache and reaction-parallel
Step-2 scoring remain unchanged.

The CORDA-like model cache is written under
`model_cache/meta_module_gem/corda_like` to prevent accidental reuse of a
FASTCORE cache.

## Audit fields

Each CORDA-like union GEM records:

- reaction-level evidence score and class;
- RNA and multiome percentile ranks;
- regulatory-support fraction;
- OT/NC support penalty;
- counts of HC, module MC, evidence MC, OT and NC reactions;
- counts of OT and NC reactions actually used as support;
- initial and final counts of parent-feasible HC/MC directions;
- completion iterations and directional closure diagnostics.

The Layer 2 result reports `params$model_completion` and a
`completion_contract`. The model-cache summary reports the evidence and support
counts for every cell-type-by-medium model.
