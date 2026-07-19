# Biological and mathematical correctness contract

This document now follows the RegCompassR 1.7.0 canonical architecture. Earlier
strict `condition × sample × cell type` metacell grouping and the independent
Pando reaction-confidence penalty are retired from the main workflow.

## Canonical architecture

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells
→ condition × cell type Pando GRNs
→ condition-specific meta-modules and local FASTCORE completion
→ one shared union-GEM and one shared medium
→ signed TF–ATAC regulation integrated into gene support
→ GPR reaction expression
→ one positive COMPASS-like reaction cost
→ directional two-step LP and descriptive condition comparison
```

## Gene-level multiome integration

RNA is mapped to bounded, zero-preserving support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficient signs and normalized magnitudes act on centered TF-by-ATAC
edge activity to produce \(R_{g,u}\in[-1,1]\). Regulation is integrated before
GPR aggregation:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

The transform preserves zero support, remains bounded, and applies activation
and repression symmetrically on the support log-odds scale.

## GPR and penalty

The canonical path uses no promiscuity weighting, Boltzmann minimum-biased AND
with `tau = 0.20`, and additive isozyme OR. Pando is not added as a separate
reaction-level term. The expression-derived reaction cost is

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

The cost is finite, strictly positive and monotonically decreasing in multiome
reaction expression. Structural reaction roles retain explicit fixed costs.

## Shared structural model

Condition-specific biological modules are completed locally with FASTCORE and
deduplicated into one union-GEM. Stoichiometry, reaction bounds, extracellular
medium, target set and target-flux fraction are shared across conditions.
Consequently, condition differences are attributable to evidence-derived costs
rather than different network structures.

Because biological samples are intentionally mixed before metacell
construction, pooled metacells do not have one sample identity. The canonical
output is therefore a descriptive condition-level comparison, not a
sample-level significance test.

See `v1.7.0-condition-pooled-architecture.md` for the full formula contract.
