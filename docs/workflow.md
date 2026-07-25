# RegCompassR 1.8.4 workflow

## Canonical data flow

```text
condition × cell type cells
→ Pando TF–peak–Human-GEM-gene models
→ significantly supported metabolic target genes
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

## 1. Pando targets and regulatory regions

Pando is fitted separately for each `condition × cell type` group. Its candidate target set is:

\[
T = G_{Human-GEM\ GPR} \cap G_{RNA\ assay}.
\]

For human hg38 analyses, the default `Pando::initiate_grn()` region set is:

\[
R = phastConsElements20Mammals.UCSC.hg38
\cup SCREEN.ccRE.UCSC.hg38.
\]

Both objects are loaded from the installed Pando package. A user-supplied `pando_initiate_args$regions` overrides the default. Because the bundled default is hg38-specific, non-human analyses must provide an appropriate region set.

A TF–peak–target row is retained as significant evidence when it passes the configured adjusted-P-value, absolute-estimate, and target-model-R² thresholds. Positive and negative coefficients both indicate regulatory evidence.

## 2. Supported metabolic gene sets and core reactions

For condition `c` and cell type `k`, let `E_{c,k}` be the significant Pando coefficient table. The supported Human-GEM metabolic genes are:

\[
M_{c,k} = \{g \in T : \exists (t,p,g) \in E_{c,k}\}.
\]

No shared-TF projection, target-target graph, top-k pruning, or GRN connected-component calculation is performed.

A reaction is a core reaction only when at least one complete GPR isozyme branch is represented:

\[
C_{c,k} = \{r : \exists j,\; GPR_{r,j} \subseteq M_{c,k}\}.
\]

Partially represented enzyme complexes are retained as diagnostics but are not core anchors. This group-level definition allows required subunits supported by different TFs or peaks to satisfy the same complete GPR branch.

## 3. Biological meta-modules

Annotation-defined expansion is ordered:

1. complete-GPR core reactions;
2. all reactions in each core reaction's annotated subsystem;
3. reactions sharing KEGG or Reactome reaction identifiers with the current set;
4. reactions sharing a master Rhea identifier with the current set.

No metabolite-neighbour, one-hop, currency-metabolite, or stoichiometric-adjacency expansion is performed.

The resulting biological reaction set is:

\[
B_{c,k} = C_{c,k} \cup S_{c,k} \cup D_{c,k} \cup R_{c,k}.
\]

## 4. Merged meta-module catalogue

Stage 3 deduplicates reaction IDs across all condition-by-cell-type biological meta-modules:

\[
B_{merged} = \bigcup_{c,k} B_{c,k},
\qquad
C_{merged} = \bigcup_{c,k} C_{c,k}.
\]

This operation does not apply medium constraints, test flux consistency, or run FASTCORE. It produces a **merged reaction catalogue**, not a GEM and not a union GEM.

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

## 5. Multiome reaction support

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

## 6. Medium-specific union GEM

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

## 7. Directional two-step LP

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

## 8. Structural comparison policy

Within one medium scenario:

- all conditions use the same union GEM;
- bounds and target-specific `vmax` are shared;
- condition differences arise from multiome penalties.

Across different media, global FASTCORE may select different support sets. Results therefore represent different structural contexts and should not be pooled into one ranking.

## 9. Targeted second-pass scoring

Selected core reaction IDs or genes are used only to define KEGG, Reactome, or master-Rhea mapping anchors. The mapped non-core reactions are scored only when present in every required final Stage 5 union GEM.

The second pass validates each cached model's checksum and medium identifier, then reuses its exact stoichiometric matrix and bounds. It does not rebuild a GEM and does not rerun FASTCORE.
