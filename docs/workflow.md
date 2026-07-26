# RegCompassR 1.8.8 workflow

## Canonical data flow

```text
all conditions within one cell type
→ one Pando structural TF–peak–GEM-gene candidate universe
→ condition-balanced multitask regression
→ global GRN backbone + condition-specific deviations
→ stability-selected condition sub-GRNs
→ condition-specific metabolic target genes
→ complete-GPR condition core reactions
→ one ordered subsystem/cross-reference expansion pass
→ merged biological reaction catalogue
→ RNA support modified by ATAC regulatory state
→ COMPASS-compatible GPR aggregation
→ one medium-specific union GEM
→ one global FASTCORE completion
→ directional two-step LP scoring
```

## 1. Shared Pando candidate background

For cell type `m`, Pando constructs:

\[
\mathcal U_m=\{e=(t,p,g)\},
\]

where `t` is a measured TF, `p` is a measured peak supported by a Pando regulatory region and a TF motif, and `g` is a GEM GPR target gene present in the RNA assay.

`Pando::prepare_grn_design()` is run once per cell type. Exact predictors are deduplicated by `(TF, ATAC feature, target)`. The candidate dictionary is independent of condition coefficients and pooled p-values.

Default regulatory regions are:

\[
R_{human}=phastCons\cup SCREEN\ ccRE,
\qquad
R_{mouse}=phastCons.
\]

## 2. Multitask condition coefficients

For edge `e=(t,p,g)` and cell `u`:

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target and predictor values are centred within condition, or within `condition × sample` when `sample_col` is supplied. The edge scale is shared across all conditions of the cell type.

RegCompass estimates:

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_c\delta_{e,c}=0.
\]

The symmetric condition basis is:

\[
H=I_C-\frac1C\mathbf 1\mathbf 1^T.
\]

Every condition has total observation weight one up to a common scaling:

\[
w_u\propto1/n_{c(u)}.
\]

The elastic-net fit uses stronger default shrinkage for condition deviations than for the global backbone. `alpha` must be below one so the centred redundant deviation representation has a unique ridge-regularised solution.

## 3. Stability-selected condition sub-GRNs

Repeated stratified subsampling gives selection frequency `Pi` and conditional sign stability `rho`:

\[
\widetilde\theta_{e,c}
=
\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

The active edge table retains global, deviation, effective and stable coefficients. Multitask coefficients do not receive classical adjusted p-values; `padj` is `NA` and the evidence type records stability selection.

The condition metabolic target set is:

\[
G_{m,c}=\{g:\exists e=(t,p,g)\text{ active in }c\}.
\]

Positive and negative active edges both establish regulatory membership.

## 4. Complete-GPR core reactions

For reaction `r` with alternative GPR branches `B_rk`:

\[
Core_{r,m,c}=1
\iff
\exists k:B_{r,k}\subseteq G_{m,c}.
\]

Required subunits joined by AND must all be present. Alternative isozyme branches remain OR alternatives. Partially represented complexes are diagnostic only.

## 5. Biological meta-modules

For each condition and cell type, annotation expansion is executed exactly once:

1. complete-GPR core reactions;
2. all reactions in core-reaction subsystems;
3. direct KEGG or Reactome reaction equivalents;
4. direct master-Rhea equivalents;
5. stop.

The resulting condition sets are merged by reaction ID:

\[
B_{merged}=\bigcup_{m,c}B_{m,c},
\qquad
C_{merged}=\bigcup_{m,c}Core_{m,c}.
\]

This is a reaction catalogue, not a GEM. Stage 3 applies no medium and runs no FASTCORE.

## 6. ATAC regulatory projection and RNA support

For each fitted edge, the ATAC-only projection coefficient is proportional to:

\[
\omega_{tpg,c}
=
\widetilde\theta_{tpg,c}\bar T_{t,m}/s_{tpg,m}.
\]

TFs sharing the same measured peak and target are signed-summed. A single target denominator is shared by all conditions:

\[
L_{g,m}=\max_c\sum_p|\psi_{p,g,c}|.
\]

The bounded modifier is:

\[
R_{g,u,c}
=
q_{g,m}
\operatorname{clip}
\left(
\frac{\sum_p\psi_{p,g,c}D_{p,u}}
{L_{g,m}+\varepsilon},-1,1
\right),
\]

where `D` is the robust ATAC deviation and

\[
q_{g,m}=\sqrt{\operatorname{clip}(R^2_{CV},0,1)}.
\]

RNA logCPM is converted to bounded support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

The regulatory modifier acts on the RNA-support log odds:

\[
C^{MO}_{g,u,c}=
\frac{C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}
{1-C^{RNA}_{g,u,c}+C^{RNA}_{g,u,c}2^{\alpha R_{g,u,c}}}.
\]

When no active edge exists for a condition, `R = 0` and `C_MO = C_RNA` exactly.

## 7. GPR reaction support and penalty

For one GPR AND branch, the canonical aggregation is:

\[
Q_{r,k,u}=\min_{g\in B_{r,k}}C^{MO}_{g,u}.
\]

`median` and `mean` are available as sensitivity options. Isozyme OR branches are additive:

\[
E_{r,u}=\sum_kQ_{r,k,u}.
\]

Reaction support becomes the LP cost:

\[
p_{r,u}=\frac1{1+\log_2(1+E_{r,u})}.
\]

## 8. Shared medium-specific union GEM

For each medium `q`, Stage 5 starts from `B_merged`, applies the medium to the parent GEM, and performs one global add-only FASTCORE completion:

\[
F_q=FASTCORE(P_q,B_{merged},C_{merged}),
\]

\[
U_q=B_{merged}\cup F_q.
\]

Only `U_q` is called a union GEM. For every condition in medium `q`:

\[
S_c=S^{U_q},
\qquad
lb_c=lb^{U_q},
\qquad
ub_c=ub^{U_q}.
\]

Therefore condition differences arise from evidence-derived penalties, not different stoichiometric models.

## 9. Directional two-step LP

For target reaction `r` and direction `d`, Step 1 computes:

\[
v^{max}_{r,d}=\max_v dv_r
\]

subject to:

\[
Sv=0,
\qquad
lb\le v\le ub.
\]

Step 2 minimizes evidence-weighted absolute flux:

\[
P^*_{r,d,u}=\min_{v,z}\sum_jp_{j,u}z_j
\]

subject to:

\[
Sv=0,
\quad
-z_j\le v_j\le z_j,
\quad
dv_r\ge\omega v^{max}_{r,d}.
\]

The default is `omega = 0.95`. Lower normalized penalty indicates stronger support in the fixed union-GEM context.

## 10. Interpretation boundary

The GRN and RNA support are learned from the same paired multiome dataset. Condition/sample centring and ATAC-only downstream projection reduce direct duplicate weighting but do not create statistically independent evidence. Condition-level metacell scores remain descriptive pseudo-observations unless biological-replicate-level validation is performed.

See [multitask GRN mathematics and object contracts](multitask-shared-grn.md) for the detailed Stage 1 derivation and table schemas.
