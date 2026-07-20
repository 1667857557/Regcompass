# Biological and mathematical correctness contract

This document follows the RegCompassR 1.7.0 canonical architecture. Earlier
strict `condition × sample × cell type` metacell grouping, Q95 calibration,
independent Pando reaction-confidence penalties and metabolite-neighbour
expansion are retired from the main workflow.

## Canonical architecture

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells with sample composition retained
→ condition × cell type Pando GRNs
→ complete-GPR core reactions
→ core-reaction subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ one shared union-GEM and one shared medium
→ coefficient-weighted ATAC accessibility integrated into RNA gene support
→ GPR reaction expression
→ one positive COMPASS-like reaction cost
→ directional two-step LP and descriptive condition comparison
```

## Gene-level multiome integration

RNA is mapped to bounded, zero-preserving support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficients are fitted from RNA+ATAC at the condition-by-cell-type level.
At per-metacell scoring, coefficient signs and normalized magnitudes act only on
robustly standardized peak-accessibility deviations to produce
\(R^{ATAC}_{g,u}\in[-1,1]\). Metacell TF RNA is not multiplied into this state.
Regulation is integrated before GPR aggregation:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}.
\]

The transform preserves zero support, remains bounded, and applies activation
and repression symmetrically on the support log-odds scale. Coefficients fitted
from the same pooled dataset remain learned parameters rather than independent
validation evidence.

## GPR and biological meta-modules

The canonical path uses no promiscuity weighting, Boltzmann minimum-biased AND
with `tau = 0.20`, and additive isozyme OR. A core reaction requires at least one
complete GPR isozyme group.

Annotation-defined biological membership may then include only:

1. reactions in each core reaction's subsystem;
2. reactions sharing KEGG or Reactome reaction identifiers;
3. reactions sharing the same master Rhea identifier.

There is no metabolite-neighbour, one-hop or generic stoichiometric-adjacency
expansion. Local FASTCORE is the only procedure that may add non-annotated
reactions, and those additions remain separately classified as feasibility
support.

## Penalty

Pando is not added as a separate reaction-level term. The expression-derived
reaction cost is

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

The cost is finite, strictly positive and monotonically decreasing in multiome
reaction expression. Structural reaction roles retain explicit fixed costs.

## Shared structural model

Condition-specific annotation-defined modules are completed locally with
FASTCORE and deduplicated into one union-GEM. Stoichiometry, reaction bounds,
extracellular medium, target set and target-flux fraction are shared across
conditions. Consequently, condition differences are attributable to
evidence-derived costs rather than different network structures.

Because biological samples are intentionally mixed before metacell
construction, pooled metacells do not have one sample identity. The canonical
output is therefore a descriptive condition-level comparison, not a
sample-level significance test.

See `v1.7.0-condition-pooled-architecture.md` for the full formula contract.
