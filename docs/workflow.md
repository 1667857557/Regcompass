# RegCompassR 1.8.3 workflow

## Canonical data flow

```text
single-cell RNA+ATAC object
→ Pando GRN per condition × cell type
→ label-guided SuperCell2 metacells within condition
→ complete-GPR core reactions
→ subsystem + KEGG/Reactome + master-Rhea meta-module expansion
→ local FASTCORE feasibility completion
→ one shared global union GEM or the validated full GEM
→ RNA support modified by accessibility-weighted Pando coefficients
→ GPR reaction expression
→ directional COMPASS-like minimum-penalty LP
→ reaction ranking, annotation and descriptive condition comparison
```

Condition is the only hard metacell stratum. `celltype_col` guides SuperCell2 before aggregation and is audited afterwards from member-cell composition. Biological-sample metadata are retained as provenance but are not used for balancing or weighting.

## GRNs and meta-modules

Pando is fitted separately for each condition × cell-type group. Significant TF–peak–gene coefficients provide signed regulatory weights estimated from the paired RNA+ATAC dataset; they are model parameters rather than independent validation evidence.

GRNs are projected to metabolic genes. A reaction is core only when at least one complete GPR isozyme group is covered. The Stage 3 biological meta-module includes:

- core-reaction subsystem members;
- reactions sharing KEGG or Reactome reaction IDs;
- reactions sharing master-Rhea IDs.

No metabolite-neighbour expansion is used. Local FASTCORE adds only reactions required for flux feasibility and labels them separately from biological members.

Completed modules are deduplicated into one global union. `model_mode = "meta_module_gem"` scores targets in this shared union; `model_mode = "full_gem"` uses the complete validated GEM.

## Multiome reaction evidence

RNA support is

\[
C^{RNA}_{g,u}=x_{g,u}/(x_{g,u}+h).
\]

Signed Pando coefficients weight standardized peak-accessibility deviations. Metacell TF RNA is not multiplied into this modifier. The bounded state is applied on the support log-odds scale:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R^{ATAC}_{g,u}}}.
\]

Protein complexes use normalized Boltzmann soft-min AND with `tau = 0.20`; isozyme groups are summed. Reaction expression becomes the positive LP cost

\[
p_{r,u}=1/[1+\log_2(1+E^{MO}_{r,u})].
\]

Only exchange, demand, sink and artificial-support reactions receive fixed structural costs.

## LP and ranking

For each target direction, step 1 maximizes directional target flux under `S v = 0` and model bounds. Step 2 requires at least `omega × vmax` and minimizes evidence-weighted absolute network flux. The structural model is shared across metacells; only reaction penalties change.

Raw minimum penalty is the primary same-target comparison. Cross-reaction priority uses `penalty / (omega × vmax)`. The workflow verifies target `vmax` invariance across metacells before descriptive condition comparison.

## Optional direct database-linked scoring

After the original core Layer 2 run, `rc_regcompass_step_target_union()` can use selected cores as mapping anchors. It directly identifies reactions sharing a KEGG, Reactome, or master-Rhea ID with each anchor and scores only linked reactions that were not original global core targets.

This optional step does **not** use same-subsystem expansion, recursive propagation, metabolite neighbours, or generic union membership. It reuses the exact cached union-GEM files and bounds from the original Layer 2 run.

## Main outputs

- `metacells`: aggregated counts, membership and composition diagnostics;
- `layer1`: RNA support, ATAC modifier, GPR diagnostics and reaction expression;
- `global_grn_meta_modules`: core, biological expansion and FASTCORE support membership;
- `microcompass`: directional `vmax`, feasibility and raw minimum penalties;
- `reaction_ranking`: within-condition reaction priorities;
- `reaction_catalog` and `reaction_evidence`: biological annotation and evidence provenance;
- optional target-union outputs: anchor-to-reaction database mappings and direct non-core LP scores.
