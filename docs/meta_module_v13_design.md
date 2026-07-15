# RegCompassR 1.4 global-metacell architecture

## Scope

The supported integrated workflow uses two explicit phases:

1. strict-stratum multiome, Pando and meta-module inference;
2. global expression-capacity recalibration, one shared structural GEM and per-metacell directional LP scoring.

RegCompassR exposes two Layer 2 structural modes:

1. `meta_module_gem`: the union of all retained strict-stratum Pando meta-modules, completed once per medium scenario with add-only FASTCORE support reactions;
2. `full_gem`: the complete medium-constrained Human-GEM, shared by all metacells.

Sample-specific or condition-specific structural models are not part of the integrated workflow because they create different feasible spaces and prevent direct score comparison.

## Phase 1: one worker per strict biological stratum

A strict stratum is defined as:

```text
condition × sample × cell type
```

Each retained stratum is processed by one upstream worker. That worker performs, in order:

```text
single-cell RNA+ATAC subset
→ SuperCell2 metacell construction
→ fragment aggregation to the metacell level
→ Pando internal peak-gene and TF-gene modeling
→ Pando GRN inference
→ Pando-derived metacell reaction confidence
→ GRN-to-GPR reaction mapping
→ subsystem / KEGG / Reactome / master-Rhea expansion
→ stratum artifact
```

Pando is the sole peak-gene model in the integrated workflow. RegCompassR does not run a second Signac `LinkPeaks()` pass.

### Pando-derived ATAC confidence

For every significant Pando target-region coefficient, RegCompassR matches the Pando region to the normalized ATAC peak matrix saved by Pando. Let `a(p,u)` be TF-IDF accessibility for peak `p` in metacell `u`. Each peak is transformed across the metacells of its strict stratum with the same robust sigmoid transform used for gene evidence:

```text
z(p,u) = [a(p,u) - median_u a(p,u)] / robust_scale_p
s(p,u) = sigmoid(z(p,u))
```

A peak with zero or non-finite accessibility receives score zero. For target gene `g`, the Pando confidence is a weighted average over matched significant regions:

```text
F_gene(g,u) = sum_p w(g,p) s(p,u) / sum_p w(g,p)
```

The default link weight is based on the absolute Pando coefficient and, when available, model fit quality. Activating and repressing coefficients both provide evidence that a regulatory region is associated with the target; coefficient sign is therefore not used as a confidence direction.

Genes without a matched significant Pando region are treated as missing ATAC confidence, not as zero-confidence observed genes. Gene confidence is aggregated through the Human-GEM GPR structure:

- AND subunits use a soft minimum;
- OR isoenzymes use the best supported complete group;
- incomplete required complexes do not define supported reaction confidence.

The resulting matrix is:

```text
reaction × metacell Pando confidence
```

It modifies only the Layer 2 penalty. It does not alter stoichiometry, reaction bounds or medium constraints.

## All-strata barrier

The global phase is blocked until all retained strict strata have completed successfully. The barrier verifies:

- every retained `condition × sample × cell-type` stratum returned an artifact;
- every artifact contains RNA metacell expression, Pando outputs, reaction confidence and meta-module membership;
- every biological sample remains represented;
- artifact identifiers match the required stratum order.

If any retained stratum fails, the upstream worker pool is released, a barrier diagnostic is written and the workflow stops. It never continues with a successful subset.

## Phase 2A: global expression-capacity recalibration

Expression capacity is not taken from stratum-local normalized outputs. After the barrier, RegCompassR combines Human-GEM GPR-gene logCPM across all metacells and recomputes the entire expression-capacity chain:

```text
all-metacell GPR-gene logCPM
→ one robust gene-score scale across all metacells
→ GPR AND/OR reaction capacity
→ one reaction-wise Q95 scale across all metacells
```

For gene `g` and metacell `u`:

```text
E(g,u) = sigmoid([x(g,u) - median_u x(g,u)] / robust_scale_g)
```

GPR reaction capacity uses the configured bottleneck-aware AND rule and isoenzyme-aware OR rule. The final relative capacity is:

```text
C_rel(r,u) = min(1, C_raw(r,u) / [Q95_r + epsilon])
```

The same reaction-specific denominator is used for every sample, condition, cell type and metacell.

## Phase 2B: global meta-module union

For strict stratum `s` and GRN module `g`, let `C(s,g)` be hard-core reactions supported by complete GPR groups. The biological envelope `B(s,g)` is obtained by ordered expansion:

```text
C(s,g)
→ reactions in core-reaction subsystems
→ reactions connected through KEGG/Reactome identifiers
→ reactions connected through master-Rhea identifiers
```

The global sets are:

```text
C_global = union over all retained s,g of C(s,g)
B_global = union over all retained s,g of B(s,g)
```

All source tables retain their original `condition`, `sample`, `cell type`, `group_id` and module provenance. The structural model uses one canonical identity:

```text
sample_id = global
module_id = GLOBAL_UNION
```

## Shared parent model and medium

For each medium scenario:

1. apply the medium constraints to Human-GEM;
2. disable demand, sink and artificial-support reactions during reconstruction;
3. verify steady-state feasibility;
4. run directional FASTCC consistency screening;
5. fix inconsistent reactions to zero bounds.

Condition-specific medium rows are rejected by the integrated workflow. A condition-specific medium would create different bounds and therefore different feasible spaces.

## Add-only FASTCORE completion

The completed global model is:

```text
M(global, medium) = B_global union support(global, medium)
```

FASTCORE support reactions are selected outside the biological envelope. The biological envelope is retained add-only. Reverse requests are handled by orienting the signed stoichiometric column and bounds rather than splitting a reaction into artificial forward and reverse copies.

A target direction can have four outcomes:

| Status | Meaning |
|---|---|
| `already_feasible` | Feasible in the global biological envelope. |
| `fastcore_completed` | Parent-feasible and restored by added support reactions. |
| `parent_blocked` | Direction allowed by bounds but infeasible in the consistent parent GEM. |
| `unresolved` | Parent-feasible but still infeasible after completion. |

With `strict = TRUE`, unresolved parent-feasible targets stop model construction. Parent-blocked targets remain explicit diagnostics and are not repaired with artificial boundaries.

## Shared-model directional microCOMPASS

Every metacell uses the same `S`, `lb`, `ub`, target set and medium-specific model. Biological differences enter through a metacell-specific penalty vector:

```text
p(r,u) =
  w_expr [-log C_rel(r,u)]
  + w_conf [-log F_pando(r,u)]
  + w_missing I_missing(r,u) M
  + optional GPR-missing term
```

For target reaction `r` and direction sign `d`, the first LP computes:

```text
maximize d × v_r
subject to S v = 0
           lb <= v <= ub
```

The second LP computes the minimum-penalty flux state while preserving a fraction `omega` of the target maximum:

```text
minimize sum_i p(i,u) a_i
subject to S v = 0
           lb <= v <= ub
           -a_i <= v_i <= a_i
           a_i >= 0
           d × v_r >= omega × vmax
```

The Layer 2 task unit is:

```text
shared model / medium × metacell
```

A worker loads the shared model once for one metacell and evaluates all target directions with that metacell's penalty vector.

## Parallel lifecycle

The workflow uses two separate parallel passes:

1. strict-stratum upstream workers;
2. shared-model-by-metacell LP workers.

After the upstream barrier, the upstream `BiocParallel` backend is explicitly stopped and large temporary objects are removed before the global stage. A fresh backend is created for Layer 2. Cleanup is executed on both success and error paths.

## Interpretation limits

- `Sv=0` is a steady-state feasibility assumption, not evidence of measured flux.
- FASTCORE is not an exact minimum-cardinality MILP.
- Pando association does not establish causal regulation.
- Reaction capacity and directional microCOMPASS scores are relative multiome-supported potentials.
- Metacells are computational units, not independent biological replicates; differential inference must operate at the biological-sample level.
- Results remain conditional on Human-GEM annotations, GPR mapping, Pando model quality, medium definition, reaction roles and numerical tolerances.
