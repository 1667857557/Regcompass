# Biological and mathematical correctness contract

This document follows the RegCompassR 1.7.0 condition-pooled architecture.
Earlier strict `condition × sample × cell type` metacell grouping, Q95
calibration, independent Pando reaction-confidence penalties, sample-level
inference APIs, and metabolite-neighbour expansion are retired from the main
workflow.

## Canonical evidence architecture

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells with sample composition retained
→ condition × cell type Pando GRNs
→ complete-GPR core reactions
→ core-reaction subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ shared union meta-module GEM or shared full GEM
→ coefficient-weighted ATAC accessibility integrated into RNA gene support
→ GPR reaction expression
→ one positive COMPASS-like reaction cost
→ directional two-step LP
→ reaction ranking and optional descriptive condition comparison
```

The workflow accepts one or more samples per condition and one or more conditions.
A single condition returns a metabolic reaction ranking without a between-group
contrast.

## Gene-level multiome integration

RNA is mapped to bounded, zero-preserving support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficients are fitted from RNA+ATAC at the condition-by-cell-type level.
At per-metacell scoring, coefficient signs and normalized magnitudes act on
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

The canonical path uses no promiscuity weighting, normalized Boltzmann soft-min
AND with `tau = 0.20`, and additive isozyme OR. A core reaction requires at least
one complete GPR isozyme group.

Annotation-defined biological membership may then include only:

1. reactions in each core reaction's subsystem;
2. reactions sharing KEGG or Reactome reaction identifiers;
3. reactions sharing the same master Rhea identifier.

There is no metabolite-neighbour, one-hop or generic stoichiometric-adjacency
expansion. Local FASTCORE is the only procedure that may add non-annotated
reactions, and those additions remain separately classified as feasibility
support.

## Penalty

Pando is not added as a separate reaction-level term. For expression-linked
reactions, the reaction cost is

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

The cost is finite, strictly positive and monotonically decreasing in multiome
reaction expression. Fixed structural costs are restricted to exchange, demand,
sink, and artificial-support reactions. Transport and cofactor reactions with
expression evidence remain governed by the integrated expression cost.

## Shared structural model modes

`meta_module_gem` deduplicates condition-specific locally completed modules into
one shared union GEM. `full_gem` retains the complete validated GEM. Both modes
share reaction penalties, medium, target set, bounds, and target-flux fraction
across all metacells.

Consequently, condition differences within one selected mode are attributable to
evidence-derived costs rather than condition-specific network structures. The two
modes are alternative structural contexts and should be reported separately.

## Ranking and inference semantics

Because biological samples are intentionally mixed before metacell construction,
pooled metacells do not retain one sample identity. The canonical outputs are:

- reaction priority ranks within every condition and cell type;
- descriptive pairwise priority differences when multiple conditions exist.

They are not biological-sample-level significance tests. With one condition, the
workflow remains valid and returns only the within-condition ranking.

See `v1.7.0-condition-pooled-architecture.md` for the full formula contract.
