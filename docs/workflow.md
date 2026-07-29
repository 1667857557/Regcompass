# RegCompassR 2.1.0 workflow

## Canonical data flow

```text
paired RNA+ATAC cells
→ condition-aware GRN inference within each broad cell type
→ condition × broad-cell-type metacells
→ supported metabolic genes and complete-GPR core reactions
→ merged biological reaction catalogue
→ regulatory modification of metacell RNA support
→ GPR-based reaction penalties
→ one shared medium-specific metabolic model
→ directional COMPASS-like LP scoring
→ reaction ranking and condition comparison
```

The workflow separates two roles:

- **Pando** estimates condition-aware regulatory models and outer-heldout target-gene scores.
- **RegCompass** maps those results to metabolic genes, reactions, metacells, penalties, and flux-constrained reaction scores.

The main comparison unit is the same reaction, direction, medium, and broad cell
type across conditions. Reference-condition coefficient contrasts are retained
for network interpretation but are not the primary Stage 4 penalty input.

## Stage 1: condition-aware GRN inference

`rc_regcompass_step_grn()` calls Pando once per broad cell type using cells from
all eligible conditions in that type. The target set is

\[
T = G_{GEM\ GPR} \cap G_{RNA\ assay}.
\]

Pando uses one shared candidate-edge dictionary and an equal-condition
coordinate system. Sparse supports are selected jointly across conditions,
while coefficients may differ in magnitude, sign, or activity status. Nested
condition-stratified outer folds generate held-out target-gene projections;
the full-data refit is retained for network interpretation only.

The canonical `candidate_screen = "motif_domain"` mode keeps motif/domain
candidates without response-dependent marginal screening. The optional
`pooled_within_condition` mode is a sensitivity analysis and is ineligible for
primary penalty construction.

Stage 1 retains two distinct coefficient views:

- absolute condition coefficients `beta_condition`, used to define supported metabolic targets and to generate the primary OOF projection;
- reference contrasts `contrast = beta_condition - beta_reference`, used only for interpretation tables.

An active coefficient can be positive or negative. Pando does not emit adjusted
P values for these regularized coefficients.

Key outputs are:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$condition_fit_status
```

## Stage 2: condition-stratified metacells

`rc_regcompass_step_metacells()` constructs SuperCell metacells separately
within each condition × broad-cell-type stratum. A metacell is rejected if its
membership mixes either condition or broad cell type.

The cache contract records the ordered cells, assay contents, reduction
embeddings, dimensions, seed, gamma, and construction thresholds. RNA and ATAC
counts are aggregated from the same membership.

The primary Pando projection is computed at the single-cell level first and is
then averaged over the exact SuperCell membership. RegCompass does not rebuild
the regulatory predictor from metacell means.

## Stage 3: metabolic targets and reaction catalogue

`rc_regcompass_step_meta_modules()` maps active absolute condition coefficients
to supported metabolic target genes. For condition `c` and cell type `k`, let
`E_{c,k}` be the active coefficient table. The supported metabolic genes are

\[
M_{c,k} = \{g \in T : \exists e \in E_{c,k},\ target(e)=g\}.
\]

A reaction is a core reaction only when at least one complete GPR isozyme branch
is represented:

\[
C_{c,k} = \{r : \exists j,\ GPR_{r,j} \subseteq M_{c,k}\}.
\]

Partially represented enzyme complexes remain diagnostics and are not core
anchors.

Each core set is expanded once in this order:

1. reactions in the same annotated subsystem;
2. direct KEGG or Reactome reaction equivalents;
3. direct master-Rhea reaction equivalents.

No newly added reaction restarts an earlier expansion step. The condition- and
cell-type-specific sets are then merged into one deduplicated reaction
catalogue. This catalogue is not yet a flux-consistent GEM.

Key outputs are:

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

## Stage 4: regulatory modification and reaction penalties

`rc_regcompass_step_layer1()` combines metacell RNA support with the Pando
outer-heldout common-support projection.

### RNA support

Metacell RNA counts are converted to latent log expression and then to bounded
gene support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Observed structural zeros remain zero; unavailable genes remain unavailable.

### Pando regulatory score

For held-out cell `i`, condition `c`, and target gene `g`, Pando stores the
absolute-condition common-support projection

\[
G^{OOF}_{i,g,c}=\sum_{e:\,target(e)=g,\,m_e=1}
\beta^{(-k)}_{e,c}z^{(-k)}_{i,e},
\]

where `m_e` is the pairwise-common or global-common estimability mask. The cell
scores are averaged within each metacell:

\[
G^{OOF}_{u,g,c}=\frac{1}{|\mathcal M_u|}
\sum_{i\in\mathcal M_u}G^{OOF}_{i,g,c}.
\]

The primary route does not use the reference contrast `Delta beta`. A
condition-estimable projection is retained only as a diagnostic comparator.

### Reliability and calibration

Pando supplies one pooled outer-heldout target `R^2`:

\[
q_g=\sqrt{clamp(R^2_{OOF,pooled,g},0,1)}.
\]

RegCompass computes one robust projection scale for each target gene and broad
cell type using all finite common-support metacell projections in that type:

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},
10^{-6}
\right).
\]

The signed modifier is

\[
R_{g,c,u}=q_g\tanh\left(\frac{G^{OOF}_{u,g,c}}{\sigma_{g,t}}\right).
\]

This is a pooled broad-cell-type calibration, not a condition-specific
rescaling. Because the scale is estimated from the projection distribution,
the modifier is not 1-homogeneous in the Pando coefficients.

### Integrated gene and reaction support

The regulatory modifier acts on RNA support on the log-odds scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

Regulatory evidence can raise or lower existing RNA support but cannot create
support when RNA support is exactly zero.

For a GPR AND group `A_{r,j}`, the default aggregation is the limiting required
subunit:

\[
Q_{r,j,u}=\min_{g\in A_{r,j}} C^{MO}_{g,u}.
\]

`"median"` and `"mean"` are available as sensitivity options. Isozyme OR
branches are additive:

\[
E_{r,u}=\sum_j Q_{r,j,u}.
\]

Reaction support is converted to the LP penalty

\[
p_{r,u}=\frac{1}{1+\log_2(1+E_{r,u})}.
\]

Missing reaction support remains unavailable rather than being converted to
zero.

## Stage 5: shared metabolic model and directional LP

`rc_regcompass_step_layer2()` applies each medium scenario to the same merged
reaction catalogue. A single global FASTCORE completion creates one
medium-specific flux-consistent model, which is reused for every condition and
metacell in that medium.

FASTCORE uses structural feasibility only; it does not use RNA or ATAC evidence.
It may add exchange, transport, cofactor-regeneration, redox-balancing, or other
support reactions required for core feasibility.

For target reaction `r` and direction `d`, the first LP computes the maximum
feasible target flux:

\[
v^{max}_{r,d}=\max_v d v_r
\]

subject to

\[
Sv=0,\qquad l\le v\le u.
\]

The second LP minimizes the network-wide penalty while requiring the target to
reach a fraction `omega` of that maximum:

\[
P^*_{r,d,u}=\min_{v,t}\sum_j p_{j,u}t_j
\]

subject to

\[
Sv=0,\quad -t_j\le v_j\le t_j,\quad
 d v_r\ge \omega v^{max}_{r,d}.
\]

The default is `omega = 0.95`. Cross-reaction ranking uses

\[
\widetilde P_{r,d,u}=\frac{P^*_{r,d,u}}
{\omega v^{max}_{r,d}}.
\]

Lower normalized penalty indicates stronger support within the fixed shared
model and medium. The score is a network-constrained support cost, not a direct
measurement of reaction flux.

## Stage 6: results and condition comparisons

`rc_regcompass_step_results()` attaches reaction annotations, evidence records,
rankings, and condition contrasts.

Within a medium scenario, all conditions use the same reaction set,
stoichiometric matrix, bounds, target direction, and target-specific `vmax`.
Condition differences therefore enter through the metacell penalty vectors.
Different media may produce different completed models and should be treated as
separate structural contexts.

`rc_test_condition_reactions()` compares the same reaction, direction, medium,
and broad cell type across conditions using metacells as statistical units. The
reported Wilcoxon and Kruskal–Wallis P values quantify within-dataset metacell
separation; they are not donor- or sample-level biological-replicate inference.

## Restartable functions

| Stage | Function | Main output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | Pando condition-aware GRNs and OOF fit contracts |
| 2 | `rc_regcompass_step_metacells()` | pure condition × broad-cell-type SuperCells |
| 3 | `rc_regcompass_step_meta_modules()` | supported genes, core reactions, merged reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | RNA support, regulatory modifiers, reaction penalties |
| 5 | `rc_regcompass_step_layer2()` | shared medium-specific models and directional LP scores |
| 6 | `rc_regcompass_step_results()` | annotations, rankings, provenance, and contrasts |
