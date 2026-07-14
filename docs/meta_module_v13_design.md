# RegCompassR 1.3 architecture and mathematical specification

## Scope

RegCompassR exposes exactly two Layer 2 structural models:

1. `full_gem`: the complete medium-constrained Human-GEM.
2. `meta_module_gem`: a sample-specific biological reaction envelope plus add-only FASTCORE support reactions.

The following v1.2 concepts are removed from the public architecture:

- target k-hop micro-GEMs;
- static module-meso-GEMs;
- automatic fallback between local model builders;
- optional support-stability enumeration;
- alternative reference-vmax modes.

## Biological reaction envelope

For sample `s` and GRN module `g`, let `C(s,g)` be direct reactions obtained by mapping Pando metabolic genes to Human-GEM GPR genes.

The biological reaction set `B(s,g)` is constructed by ordered expansion:

```text
C(s,g)
→ all reactions in core-reaction subsystems
→ all reactions in subsystems linked by shared KEGG/Reactome reaction IDs
→ all reactions in subsystems linked by shared master-Rhea IDs
```

`B(s,g)` is retained in the final local model and is never penalized by FASTCORE LP-10. The mandatory directional target set remains `C(s,g)`. Peripheral members of `B(s,g)` are biological context, not mandatory hard-core reactions.

This distinction prevents a blocked peripheral annotation from invalidating an otherwise usable GRN-defined model while preserving every biologically selected reaction for downstream inspection.

## Parent model

For each medium scenario `m`:

1. Apply medium constraints to the complete Human-GEM.
2. Disable reactions annotated as `demand`, `sink`, or `artificial_support`.
3. Solve a zero-objective steady-state LP to verify that the bounded model is feasible.
4. Run directional FASTCC screening.
5. Fix inconsistent reactions to zero bounds.

The parent model is shared by all sample-module tasks under the same medium scenario.

No artificial exchange or gap-filling reaction is introduced during completion.

## Signed reaction directions

Human-GEM uses one signed flux variable per reaction:

```text
lb < 0 < ub     reversible
0 <= lb < ub    forward-only
lb < ub <= 0    reverse-only
lb = ub = 0     blocked
```

Positive flux follows the stored stoichiometric column. Negative flux follows the reverse direction.

A reverse core task is temporarily oriented as a positive task:

```text
S_j'  = -S_j
lb_j' = -ub_j
ub_j' = -lb_j
```

The operation is a bijective variable transformation. It avoids splitting a reversible reaction into separate positive and negative variables, which could create a false internal cycle.

## FASTCC screening

For a candidate set `J`, the LP-7 form is:

```text
maximize sum_j z_j
subject to S v = 0
           lb <= v <= ub
           v_j - z_j >= 0
           0 <= z_j <= epsilon
```

Reactions with `z_j >= epsilon` are consistent in the tested orientation. Remaining reversible candidates are tested after stoichiometric and bound orientation. Repeated LP-7 calls plus singleton checks classify all parent reactions.

`fastcore_epsilon` therefore defines both the numerical consistency threshold and the minimum core flux used during reconstruction.

## Add-only FASTCORE completion

For unresolved parent-feasible core targets `J` and the currently penalized reaction set `P`, RegCompass runs LP-7 followed by LP-10.

### LP-7

LP-7 identifies a subset `K` of unresolved core reactions that can be simultaneously active:

```text
maximize sum_j z_j
subject to S v = 0
           lb <= v <= ub
           v_j - z_j >= 0
           0 <= z_j <= epsilon
```

### LP-10

LP-10 minimizes absolute flux through reactions outside the biological envelope and outside previously selected support:

```text
minimize sum_i a_i, i in P
subject to S v = 0
           scaled_lb <= v <= scaled_ub
           -a_i <= v_i <= a_i
           a_i >= 0
           v_k >= scaled_epsilon, k in K
```

The implementation follows the FASTCORE numerical scaling convention:

```text
scaling_factor = 1e5
scaled_epsilon = 1e5 × epsilon
scaled_lb      = 1e5 × lb
scaled_ub      = 1e5 × ub
```

A penalized reaction is selected as support when its absolute scaled-solution flux is at least the original `epsilon`.

Selected support reactions are removed from `P`, added to the current local set, and the unresolved targets are retested. Singleton LP-10 calls are used only when the batched LP-7/LP-10 step makes no progress.

The final model is:

```text
M(s,g,m) = B(s,g) union support(s,g,m)
```

This is an add-only adaptation of FASTCORE. It seeks a compact support set but does not claim an exact minimum-cardinality solution.

## Strict completion semantics

A target direction can have four outcomes:

| Status | Meaning |
|---|---|
| `already_feasible` | Feasible in the biological envelope before FASTCORE. |
| `fastcore_completed` | Parent-feasible and restored by added support reactions. |
| `parent_blocked` | Infeasible in the medium-specific consistent parent GEM. |
| `unresolved` | Parent-feasible but still infeasible after completion. |

With `strict = TRUE`, the build stops only for `unresolved` parent-feasible targets. `parent_blocked` targets remain explicit diagnostics and are never silently repaired by artificial boundaries.

## Structural provenance

Every completed model stores:

```text
reaction_meta$biological_meta_module_member
reaction_meta$fastcore_support
reaction_meta$support_only
closure_diagnostics
completion_iterations
build_params
target_status
```

The biological envelope and operational support set are therefore distinguishable in all downstream analyses.

## Cache scope

A completed meta-module model is cached once for each:

```text
sample_id × module_id × medium_scenario
```

All core reaction directions from that sample-module reuse the same completed model. A sample-specific model is scored only against Layer 2 units whose `sample_id` matches the model.

The row identifier is:

```text
sample=<sample>::module=<module>::reaction=<reaction>::direction=<direction>::medium=<medium>
```

## Directional microCOMPASS

For direction sign `d` (`+1` forward, `-1` reverse), the first LP computes:

```text
maximize d × v_target
subject to S v = 0
           lb <= v <= ub
```

The second LP uses one signed reaction variable `v_i` and one absolute-value auxiliary variable `a_i`:

```text
minimize sum_i penalty_i × a_i
subject to S v = 0
           lb <= v <= ub
           -a_i <= v_i <= a_i
           a_i >= 0
           d × v_target >= omega × vmax
```

This formulation preserves forced non-zero positive and negative bounds. It replaces the former loose `vplus - vminus` representation, in which both components could be positive and violate the intended signed bound semantics.

In `meta_module_gem` mode, `vmax` is always computed in the completed local model. No full-GEM or hybrid reference option is exposed.

## Evidence separation

RNA-GPR capacity and ATAC-supported confidence affect the microCOMPASS penalty only. They do not change FASTCORE structural membership.

FASTCORE structural selection uses:

- stoichiometry;
- reaction bounds;
- medium constraints;
- biological membership;
- core reaction directions;
- reaction-role exclusions.

This separation avoids circularly constructing a network from the same expression evidence later used to score it.

## Required validation

The package regression suite covers:

- forward-only, reverse-only, reversible, and blocked bounds;
- forced non-zero positive and negative bounds;
- preservation of every biological-envelope reaction;
- compact forward and reverse support completion;
- the FASTCORE `1e5` LP-10 scaling convention;
- parent-blocked targets without artificial gap filling;
- avoidance of false reversible split cycles;
- one cache per sample, module, and medium;
- sample-matched scoring;
- the two-mode public API contract;
- labeled row-ID parsing.

## Interpretation limits

- Flux consistency is not thermodynamic feasibility.
- FASTCORE is not an exact minimum-cardinality MILP.
- A completed model guarantees requested parent-feasible core directions, not simultaneous activation of every core direction in one flux vector.
- Biological-envelope membership does not imply that every envelope reaction must carry flux.
- microCOMPASS scores are multiome-supported reaction potentials, not measured metabolic fluxes.
- Results remain conditional on Human-GEM annotations, GPR mapping, medium definition, reaction roles, and numerical tolerances.
