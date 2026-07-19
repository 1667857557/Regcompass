# RegCompassR 1.7.0 workflow

## Canonical data flow

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 metacells
→ condition × cell type Pando GRNs
→ GRN-derived metabolic meta-modules
→ local FASTCORE completion
→ shared union-GEM
→ TF–ATAC regulation integrated into gene support
→ GPR reaction expression
→ COMPASS-like directional minimum-penalty LP
→ descriptive condition comparison
```

## Metacell scope

The canonical workflow deliberately supplies cells from every biological sample
within the same condition and cell type to one SuperCell2 run. The original
`sample_col` is required to define the input design, but pooled metacells do not
retain one sample identity. `inference_unit = "metacell"` is therefore fixed.

The v1.7.0 canonical path requires `fragment_files = FALSE` and aggregates the
existing ATAC peak-count assay. It does not silently combine sample fragment
files without an explicit per-file barcode map.

## Pando and meta-modules

Pando is run independently for every condition-by-cell-type metacell group.
Significant TF–peak–gene coefficients define regulatory direction and relative
edge weight. The same GRNs are projected to metabolic genes, mapped through
complete GPR AND groups, expanded by subsystem and database cross-references,
and completed locally with FASTCORE.

All completed condition-specific modules are deduplicated into one shared
union-GEM. The union-GEM, medium, bounds, target reactions and target-flux
fraction are identical for all conditions.

## Multiome evidence

RNA support is

\[
C^{RNA}_{g,u}=x_{g,u}/(x_{g,u}+h).
\]

The signed Pando modifier is applied on the support log-odds scale:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

The transform is bounded and zero preserving. Protein complexes use a
Boltzmann minimum-biased AND rule with `tau = 0.20`; isozyme groups are summed;
no promiscuity weighting is applied.

Reaction expression becomes one positive LP cost:

\[
p_{r,u}=1/[1+\log_2(1+E^{MO}_{r,u})].
\]

There is no independent Pando reaction-confidence penalty.

## Outputs

- `layer1`: RNA support, regulatory modifier, multiome gene support and reaction expression;
- `grn_meta_modules`: condition-specific modules and global union membership;
- `microcompass`: directional maximum flux, feasibility and minimum penalties;
- `condition_summary`: median penalty and support per condition;
- `condition_contrast`: two-condition support difference when exactly two conditions exist.
