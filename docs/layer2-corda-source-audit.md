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

Four operations are outside the original MATLAB function:

1. cell-type multiome evidence is mapped to HC, MC, NC and OT groups;
2. the requested medium is applied before reconstruction;
3. CORDA2 still performs its original directional decomposition internally, but finalization is reaction-level: if either directional copy survives, the reaction is retained and its original medium-constrained bounds are restored; all core reactions are retained unconditionally as the immutable structural backbone, even when CORDA2 marks their tested direction as blocked;
4. no post-reconstruction closure LP is run. Scoring targets are enumerated directly from the reconstructed final GEM bounds, and microCOMPASS computes the directional `vmax` exactly once as COMPASS Step 1 before the penalty LP.

The direction relaxation and hard core retention are applied only when directional CORDA2 variables are merged back to reaction space. They do not alter Step 1, Step 2.1, Step 2.2 or Step 3 dependency calculations, including the original opposite-direction closure used while assessing an individual split target.

Core structural retention and flux feasibility are deliberately separated. For every required core reaction `r`, RegCompass enforces `r` in the final reaction set. A retained core direction may still have `vmax < flux_threshold` under the final stoichiometric and medium constraints; in that case microCOMPASS records the direction as infeasible and skips the penalty-minimization LP rather than deleting the reaction.

For `target_direction = "both"`, target enumeration uses only the final reconstructed GEM. A reversible retained core reaction contributes forward and reverse candidates according to its restored final `lb/ub`; an irreversible parent reaction remains irreversible. No parent-GEM target enumeration or parent-GEM scoring LP remains after CORDA2 finalization.

## Parallelism

Within each original CORDA2 stage, RegCompass parallelizes independent directional targets with target-isolated solver state and a stage barrier before confidence-state reduction. Each worker uses one HiGHS thread. After reconstruction there is no closure worker pass; the existing microCOMPASS shared directional-`vmax` cache owns COMPASS Step 1 and reuses each result across matching metacells.

## Validation

`tests/corda-synthetic-check.R` verifies defaults, accepted and rejected parameters, reversible decomposition, opposite-direction closure, target constraint handling, Step-1 dependencies and the HC-to-MC/HC-to-NC matrices. Testthat source-contract tests additionally verify immutable core retention, restoration of reversible parent bounds when either CORDA2 split direction is selected, absence of post-reconstruction closure LP calls, final-GEM-only target enumeration, and ownership of directional feasibility by the microCOMPASS `vmax` cache.
