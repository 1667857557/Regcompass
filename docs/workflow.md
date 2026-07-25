# RegCompassR 1.8.4 workflow

## Canonical data flow

```text
condition × cell type cells
→ Pando TF–peak–metabolic-gene GRN
→ condition-level multimodal metacells
→ signed metabolic-gene projection
→ GRN connected components
→ complete-GPR core reactions
→ core-subsystem expansion
→ KEGG/Reactome reaction-equivalence expansion
→ master-Rhea reaction-equivalence expansion
→ biological meta-modules
→ merged meta-module reaction catalogue
→ RNA support modified by ATAC regulatory state
→ medium-specific union GEM
→ global FASTCORE completion
→ directional two-step COMPASS-like LP
→ within-condition ranking and descriptive contrasts
```

## 1. GRN and metacell analysis units

Pando is fitted separately for each `condition × cell type` group. Metacells are condition-level pseudo-observations with cell-type-guided construction. Original biological-sample metadata are provenance only and are not used for balancing or weighting.

## 2. Biological meta-modules

For a GRN component `m`, a reaction is a core reaction only when at least one complete GPR isozyme group is represented:

\[
C_m = \{r : \exists k,\; GPR_{r,k} \subseteq G_m\}.
\]

Partially represented enzyme complexes are retained as diagnostics but are not core anchors.

Annotation-defined expansion is ordered:

1. complete-GPR core reactions;
2. all reactions in each core reaction's annotated subsystem;
3. reactions sharing KEGG or Reactome reaction identifiers with the current set;
4. reactions sharing a master Rhea identifier with the current set.

No metabolite-neighbour, one-hop, currency-metabolite, or stoichiometric-adjacency expansion is performed.

The resulting biological reaction set is:

\[
B_m = C_m \cup S_m \cup D_m \cup R_m.
\]

## 3. Merged meta-module catalogue

Stage 3 deduplicates reaction IDs across all biological meta-modules:

\[
B_{merged} = \bigcup_m B_m,
\qquad
C_{merged} = \bigcup_m C_m.
\]

This operation does not apply medium constraints, does not test flux consistency, and does not run FASTCORE. It produces a **merged reaction catalogue**, not a GEM and not a union GEM.

```r
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

## 4. Multiome reaction support

For gene `g` in metacell `u`, RNA logCPM is converted to bounded support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Peak accessibility is robustly centred and scaled within cell type across conditions. Signed Pando coefficients define an ATAC regulatory modifier `R_{g,u}`. The modifier acts on RNA support on the log-odds scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

ATAC can raise or lower existing RNA support but cannot create support when RNA support is zero.

Complete GPR complexes use a normalized soft minimum; isozyme groups are additive. Reaction expression support is converted to the LP cost:

\[
p_{r,u}=\frac{1}{1+\log_2(1+E_{r,u})}.
\]

## 5. Medium-specific union GEM

For each medium scenario `q`, Stage 5 builds a medium-constrained FASTCC-consistent parent GEM. Demand, sink, and artificial-support reactions are disabled for reconstruction.

The initial structural set is `B_merged`, and the target set is `C_merged`. A single global add-only FASTCORE completion identifies support reactions:

\[
F_q = FASTCORE(P_q, B_{merged}, C_{merged}).
\]

The union GEM is:

\[
U_q = B_{merged} \cup F_q.
\]

Only `U_q` is called a union GEM. One cached `U_q` is shared across all conditions and metacells analysed under medium `q`.

FASTCORE does not use RNA or ATAC evidence. It adds flux-consistent support reactions required for direction-specific core feasibility under the medium. Support may include exchange, transport, cofactor-regeneration, redox-balancing, or connecting internal reactions.

Global FASTCORE controls are supplied through:

```r
layer2_args$model_params <- list(
  completion_time_limit = 600,
  fastcore_epsilon = 1e-4,
  max_support_reactions = 2000,
  strict = TRUE
)
```

## 6. Directional two-step LP

For target reaction `r` and direction `d`, Step 1 computes:

\[
v^{max}_{r,d}=\max_v d v_r
\]

subject to:

\[
Sv=0,\qquad l\le v\le u.
\]

Step 2 minimizes network-wide evidence discordance:

\[
P^*_{r,d,u}=\min_{v,t}\sum_j p_{j,u}t_j
\]

subject to:

\[
Sv=0,
\quad -t_j\le v_j\le t_j,
\quad d v_r\ge \omega v^{max}_{r,d}.
\]

The default target-flux fraction is `omega = 0.95`.

Cross-reaction ranking uses:

\[
\widetilde P_{r,d,u}=\frac{P^*_{r,d,u}}
{\omega v^{max}_{r,d}}.
\]

Lower normalized penalty means stronger support in the fixed union-GEM context.

## 7. Structural comparison policy

Within one medium scenario:

- all conditions use the same union GEM;
- bounds and target-specific `vmax` are shared;
- condition differences arise from multiome penalties.

Across different media, global FASTCORE may select different support sets. Results therefore represent different structural contexts and should not be pooled into one ranking.

## 8. Targeted second-pass scoring

Selected core reaction IDs or genes are used only to define KEGG, Reactome, or master-Rhea mapping anchors. The mapped non-core reactions are scored only when present in every required final Stage 5 union GEM.

The second pass validates each cached model's checksum and medium identifier, then reuses its exact stoichiometric matrix and bounds. It does not rebuild a GEM and does not rerun FASTCORE.
