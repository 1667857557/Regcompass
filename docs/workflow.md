# RegCompassR 1.9.1 workflow

## Canonical data flow

```text
cell type cells across conditions
→ shared Pando edge design and independent condition coefficients
→ supported metabolic target genes
→ complete-GPR core reactions
→ one ordered subsystem/cross-reference expansion pass
→ biological meta-modules
→ merged meta-module reaction catalogue
→ RNA support modified by ATAC regulatory state
→ COMPASS-compatible GPR-AND aggregation
→ medium-specific union GEM
→ global FASTCORE completion
→ directional two-step COMPASS-like LP
→ within-condition ranking and descriptive contrasts
```

## 1. Pando targets, motifs, and regulatory regions

Pando is fitted once per cell type across conditions. The complete edge
dictionary, edge eligibility mask, and pooled final-predictor scale are shared.
Within each target, the lambda path and selected lambda are also shared;
condition coefficients are independently estimated. Its candidate target set is:

\[
T = G_{GEM\ GPR} \cap G_{RNA\ assay}.
\]

When `pfm` is omitted, RegCompass loads the Pando data object:

```r
data("motifs", package = "Pando")
pfm <- motifs
```

A user-supplied motif collection overrides this default.

The default `Pando::initiate_grn()` region set depends on `species`:

\[
R_{human} = phastConsElements20Mammals.UCSC.hg38
\cup SCREEN.ccRE.UCSC.hg38,
\]

\[
R_{mouse} = phastConsElements20Mammals.UCSC.hg38.
\]

The objects are loaded from the installed Pando package. A user-supplied `pando_initiate_args$regions` overrides either species-specific default. The selected policy is recorded in `step1$grn_result$normalization_policy$pando_regions`.

An active condition TF–peak–target coefficient is retained when it passes the
configured absolute-estimate and target-model-R² thresholds. Positive and
negative coefficients both indicate regulatory evidence. The regularized
solver does not emit adjusted P values.

## 2. Supported metabolic gene sets and core reactions

For condition `c` and cell type `k`, let `E_{c,k}` be the active Pando coefficient table after the configured absolute-effect and model-quality filters. The supported GEM metabolic genes are:

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

Annotation-defined expansion is executed exactly once in the following order:

1. begin with complete-GPR core reactions `C_{c,k}`;
2. add all reactions in each core reaction's annotated subsystem;
3. from the resulting core-plus-subsystem set, add direct KEGG or Reactome reaction equivalents;
4. from the resulting core-plus-subsystem-plus-database set, add direct master-Rhea reaction equivalents;
5. stop.

A reaction introduced at the master-Rhea step does not trigger another KEGG/Reactome or subsystem pass.

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
step3$condition_modules$meta_module_summary$expansion_policy
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

## 5. Multiome reaction support and GPR aggregation

For gene `g` in metacell `u`, RNA logCPM is converted to bounded support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

For condition \(c\), RegCompass uses the explicit Pando contrast
\(\Delta\beta_{e,c}=\beta_{e,c}-\beta_{e,reference}\). It reconstructs each
edge predictor from metacell TF RNA and peak ATAC, then applies the exact pooled
center and scale stored by Pando:

\[
z_{e,u} =
\frac{\mathrm{RNA}_{TF(e),u}\mathrm{ATAC}_{peak(e),u}-\mu_e}{s_e}.
\]

The raw model-space projection and bounded modifier are

\[
M_{g,c,u} =
\sqrt{\operatorname{clamp}(R^2_{g,c},0,1)}
\sum_{e:\,target(e)=g}\Delta\beta_{e,c}z_{e,u},
\qquad
R_{g,c,u}=\tanh(M_{g,c,u}).
\]

There is no downstream metacell robust re-scaling and no division by the
absolute sum of edge effects. The modifier acts on RNA support on the log-odds
scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

ATAC can raise or lower existing RNA support but cannot create support when RNA support is zero.

For one GPR AND group `A_{r,j}`, RegCompass uses one of the three COMPASS aggregation functions:

\[
Q^{min}_{r,j,u}=\min_{g\in A_{r,j}} C^{MO}_{g,u},
\]

\[
Q^{median}_{r,j,u}=\operatorname{median}_{g\in A_{r,j}} C^{MO}_{g,u},
\]

\[
Q^{mean}_{r,j,u}=\frac{1}{|A_{r,j}|}
\sum_{g\in A_{r,j}} C^{MO}_{g,u}.
\]

The canonical default is `gpr_and_method = "min"`, representing the limiting required subunit. `"median"` and `"mean"` are available for sensitivity analysis.

Isozyme OR branches remain additive:

\[
E_{r,u}=\sum_j Q_{r,j,u}.
\]

Reaction expression support is converted to the LP cost:

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
