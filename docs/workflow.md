# RegCompassR 1.7.0 workflow

## Canonical data flow

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 metacells with biological-sample composition retained
→ one Pando GRN per condition × cell type
→ complete-GPR core reactions
→ subsystem + KEGG/Reactome + master-Rhea expansion
→ one bounded non-structural metabolite-neighbour hop
→ local FASTCORE completion
→ one shared union-GEM and shared medium
→ Pando-coefficient-weighted ATAC state integrated into RNA support
→ GPR reaction expression
→ COMPASS-like directional minimum-penalty LP
→ descriptive condition comparison
```

## Metacell scope

The canonical workflow deliberately supplies cells from every biological sample
within the same condition and cell type to one SuperCell2 run. Each sample must
map to one condition; strict defaults require at least two biological samples
per condition.

Pooled metacells do not have one biological-sample identity, but the workflow
retains cell membership and reports per-metacell sample counts, fractions,
dominant-sample fraction and effective sample number. Pooling is cell-count
weighted. `inference_unit = "metacell"` is fixed, and these units are descriptive
pseudo-observations rather than independent biological replicates.

The v1.7.0 canonical path requires `fragment_files = FALSE` and aggregates the
existing ATAC peak-count assay. It does not silently combine sample fragment
files without an explicit per-file barcode map.

## Pando and meta-modules

Pando is run independently for every condition-by-cell-type metacell group.
Significant TF–peak–gene coefficients define regulatory direction and relative
edge weight. These coefficients are fitted from RNA+ATAC and are not independent
validation evidence.

GRNs are projected to metabolic genes. A reaction is core only when at least one
complete GPR isozyme group is present. Biological expansion then adds reactions
from the core subsystem, shared KEGG/Reactome reaction IDs, shared master Rhea
IDs, and exactly one bounded non-structural metabolite-neighbour hop. High-degree
currency metabolites are excluded from that hop. Local FASTCORE subsequently
adds only reactions required for feasibility.

All completed condition-specific modules are deduplicated into one shared
union-GEM. The union-GEM, medium, bounds, target reactions and target-flux
fraction are identical for all conditions.

## Multiome evidence

RNA support is

\[
C^{RNA}_{g,u}=x_{g,u}/(x_{g,u}+h).
\]

At the metacell scoring stage, signed Pando coefficients weight robustly
standardized **peak-accessibility** deviations. Metacell TF RNA is not multiplied
into the modifier. The bounded state \(R^{ATAC}_{g,u}\in[-1,1]\) is applied on
the support log-odds scale:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}.
\]

The transform is bounded and zero preserving. Protein complexes use the
normalized, monotone Boltzmann soft-min AND rule with `tau = 0.20`; isozyme
groups are summed and no promiscuity weighting is applied.

Reaction expression becomes one positive LP cost:

\[
p_{r,u}=1/[1+\log_2(1+E^{MO}_{r,u})].
\]

There is no independent Pando reaction-confidence penalty, Q95 calibration,
confidence-alignment matrix or `penalty_weights` term.

## LP and outputs

For each target and direction, step 1 maximizes directional target flux under
`S v = 0` and model bounds. Step 2 requires at least `omega × vmax` and minimizes
the weighted absolute network flux. The stoichiometric model is shared; only the
unit-specific reaction penalty changes.

- `metacells`: pooled metacells, membership and sample composition;
- `layer1`: RNA support, ATAC-derived modifier, multiome gene support and `reaction_expression`;
- `grn_meta_modules`: biological membership, FASTCORE completion and global union membership;
- `microcompass`: directional maximum flux, feasibility and minimum penalties;
- `condition_summary` and `condition_contrast`: descriptive within-cell-type comparisons.
