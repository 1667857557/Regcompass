# Original CORDA2 source-to-R audit

Reference: `schultzdre/Constraint-Based-Modeling/CORDA2.m`.

The audit separates the RegCompass confidence adapter from the CORDA2 reconstruction state machine and downstream directional scoring.

## Source correspondence

| Original MATLAB operation | RegCompass implementation |
|---|---|
| Defaults `MCxNCthresh=2`, `constraint=1`, `constrainby='val'`, `om=1e4`, `ci=0.01` | `.rc_corda2_original_args()` |
| Split reactions satisfying `lb < 0 && ub >= 0`; suffix `_CORDA_rev_rxn` | `.rc_corda2_split_original()` |
| Set all split variables non-negative and reverse upper bound to `-lb` | `.rc_corda2_split_original()` |
| Close the reverse copy when testing forward, or close the original when testing reverse | `.rc_corda2_close_opposite()` |
| Maximize target before assigning its value or percentage constraint | `.rc_corda2_maximize_target()`, `.rc_corda2_constrain_target()` |
| Step-1 costs: other `1e-3`, MC `sqrt(om)`, NC `om` | `.rc_corda2_stage_cost(..., "stage1")` |
| Increase each newly observed high-cost dependency by `1+ci` | `.rc_corda2_dependency_assessment()` |
| Stop when the directional dependency set no longer changes | `.rc_corda2_dependency_assessment()` |
| Promote MC and NC directions required by HC | `.rc_corda_build_three_stage()` Step 1 |
| Build MC-to-NC dependency matrix | `.rc_corda_build_three_stage()` Step 2.1 |
| Promote NC directions with occurrence at least `MCxNCthresh` | `.rc_corda_build_three_stage()` Step 2.2 |
| Block remaining NC directions and test MC feasibility by maximization | `.rc_corda_build_three_stage()` Step 2.2 |
| Block remaining MC/NC and add required OT directions with cost `om` | `.rc_corda_build_three_stage()` Step 3 |
| Return `HCtoMC`, `HCtoNC`, `MCtoNC` and `rescue` | reconstruction object and union-GEM fields |
| Merge selected directional copies back to reaction identifiers | `.rc_corda2_apply_direction_bounds()` |

## Parameter boundary

The public CORDA2 parameter list contains only the five adjustable inputs present in the original function:

```r
corda2_args = list(
  MCxNCthresh = 2,
  constraint = 1,
  constrainby = "val",
  om = 1e4,
  ci = 0.01
)
```

Internal source constants such as the positive-flux threshold and Step-1 baseline cost are not public parameters.

## RegCompass adapters

Three operations are outside the original MATLAB function:

1. cell-type multiome evidence is mapped to HC, MC, NC and OT groups;
2. the requested medium is applied before reconstruction;
3. final selected directions retain the medium-constrained parent bounds instead of reopening every retained reaction to unconditional `±1000` bounds.

The third operation is required to preserve the medium supplied to Layer 2. It changes only the final bound restoration adapter, not dependency assessment or reaction inclusion.

## Parallelism

The original directional target loop is serial within each reconstruction. RegCompass parallelizes only independent `cell type × medium` models. Each worker uses one solver thread, writes its cache atomically and releases the native solver state before returning to the pool.

## Validation

`tests/corda-synthetic-check.R` verifies defaults, accepted and rejected parameters, reversible decomposition, opposite-direction closure, target constraint handling, Step-1 dependencies and the HC-to-MC/HC-to-NC matrices. Testthat source-contract tests additionally verify the direct call path and absence of function-override wrappers.
