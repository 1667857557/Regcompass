# Mathematical model

This chapter contains the mathematical definitions used by the canonical
RegCompassR workflow. Other documentation describes functions, parameters,
outputs, and interpretation without repeating these equations.

## 1. Condition-aware Pando model

For one broad cell type, target gene `g`, condition `c`, and candidate edge `e`,
the predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

### Equal-condition transform

For `K` conditions,

\[
\mu_e=\frac{1}{K}\sum_c\bar x_{e,c},
\qquad
s_e^2=\frac{1}{K}\sum_c\frac{1}{n_c}
\sum_{i\in c}(x_{e,i}-\bar x_{e,c})^2,
\]

\[
z_{e,i}=\frac{x_{e,i}-\mu_e}{s_e}.
\]

The target response uses the analogous equal-condition center and
within-condition variance scale. Outer-fold transforms are estimated from the
training cells only.

### Sparse-group multitask objective

Let `B` contain one coefficient column per condition. For fixed `lambda`,
`alpha`, and `rho = condition_mix`, Pando minimizes

\[
\sum_c\frac{1}{2n_c}\lVert y_c-a_c-X_c\beta_c\rVert_2^2
+\frac{\lambda(1-\alpha)}{2}\lVert B\rVert_F^2
+\lambda\alpha\left[(1-\rho)\sum_e\lVert B_{e,\cdot}\rVert_2
+\rho\sum_{e,c}|\beta_{e,c}|\right].
\]

The row-wise group term couples edge selection across conditions; the
entry-wise term permits condition-specific zeros. Pando then performs a
support-constrained common-metric ridge refit.

### Reference contrast

For reference condition `r`,

\[
\Delta\beta_{e,c}=\beta_{e,c}-\beta_{e,r}.
\]

This contrast is retained for interpretation. The canonical penalty uses
absolute condition coefficients.

## 2. Outer-heldout regulatory projection

For held-out cell `i`, outer fold `k`, condition `c`, target gene `g`, and common
estimability mask `m_e`,

\[
G^{OOF}_{i,g,c}=\sum_{e:\,target(e)=g,\,m_e=1}
\beta^{(-k)}_{e,c}z^{(-k)}_{i,e}.
\]

Two-condition comparisons use pairwise-common support; multi-condition
comparisons use global-common support. Condition-estimable projections are
diagnostic only.

For metacell `u` with membership set `M_u`,

\[
G^{OOF}_{u,g,c}=\frac{1}{|M_u|}\sum_{i\in M_u}G^{OOF}_{i,g,c}.
\]

The single-cell projection is computed before metacell aggregation.

## 3. Reliability and calibration

The pooled outer-heldout target fit is

\[
R^2_{OOF,g}=1-
\frac{\sum_c\sum_{i\in c}(y_{i,g}-\hat y^{OOF}_{i,g})^2}
{\sum_c\sum_{i\in c}(y_{i,g}-\bar y_{c,g})^2}.
\]

The reliability weight is

\[
q_g=\sqrt{clamp(R^2_{OOF,g},0,1)}.
\]

For target gene `g` and broad cell type `t`, RegCompass computes one pooled
robust scale over finite common-support metacell projections:

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},
10^{-6}
\right).
\]

If fewer than two finite values are available or the projection has no range,
the fallback scale is `1`. The bounded signed modifier is

\[
R_{g,c,u}=q_g\tanh\left(\frac{G^{OOF}_{u,g,c}}{\sigma_{g,t}}\right).
\]

Because the scale is estimated from the projection distribution, this modifier
is not 1-homogeneous in the Pando coefficients.

## 4. Gene support

Metacell RNA expression is converted to bounded support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

The regulatory modifier acts on the RNA-support odds:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

RNA support equal to `0` or `1` remains unchanged.

## 5. Supported genes and core reactions

For condition `c` and cell type `t`, let `E_{c,t}` be the active absolute
condition-coefficient table and `T` the metabolic target-gene set:

\[
M_{c,t}=\{g\in T:\exists e\in E_{c,t},\ target(e)=g\}.
\]

A reaction is a core reaction when at least one complete GPR isozyme branch is
represented:

\[
C_{c,t}=\{r:\exists j,\ GPR_{r,j}\subseteq M_{c,t}\}.
\]

The biological catalogue is the union of core reactions and the configured
single-pass subsystem, KEGG/Reactome, and master-Rhea expansions.

## 6. GPR aggregation and reaction penalty

For GPR AND group `A_{r,j}`, the canonical rule is

\[
Q_{r,j,u}=\min_{g\in A_{r,j}}C^{MO}_{g,u}.
\]

`median` and `mean` are optional sensitivity rules. Isozyme OR branches are
additive:

\[
E_{r,u}=\sum_jQ_{r,j,u}.
\]

Reaction support is converted to the LP penalty

\[
p_{r,u}=\frac{1}{1+\log_2(1+E_{r,u})}.
\]

## 7. Shared metabolic model

For each medium, RegCompass applies one global FASTCORE completion to the merged
reaction catalogue and reuses the resulting stoichiometric model, bounds, and
reaction order for every condition and metacell.

For target reaction `r` and direction `d`, the first LP computes

\[
v^{max}_{r,d}=\max_v d v_r
\]

subject to

\[
Sv=0,\qquad l\le v\le u.
\]

The second LP minimizes the network-wide penalty while maintaining a fraction
`omega` of the target maximum:

\[
P^*_{r,d,u}=\min_{v,t}\sum_jp_{j,u}t_j
\]

subject to

\[
Sv=0,\qquad -t_j\le v_j\le t_j,
\qquad d v_r\ge\omega v^{max}_{r,d}.
\]

The normalized penalty is

\[
\widetilde P_{r,d,u}=\frac{P^*_{r,d,u}}
{\omega v^{max}_{r,d}}.
\]

Lower normalized penalty indicates stronger network-constrained support for the
target direction; it is not a measured flux.

For a fixed model, bounds, target, direction, and `vmax`, the LP optimum is
coordinatewise nondecreasing, positively 1-homogeneous, and concave in the
reaction-penalty vector.

## 8. Condition statistics

`rc_test_condition_reactions()` transforms the normalized penalty to

\[
S_{r,d,u}=-\log(\widetilde P_{r,d,u}+\epsilon).
\]

Pairwise tests use Wilcoxon rank-sum statistics; analyses with at least three
conditions may also use a Kruskal-Wallis omnibus test. Metacells are the
statistical units, so these P values describe within-dataset metacell separation
and are not donor-level biological-replicate inference.
