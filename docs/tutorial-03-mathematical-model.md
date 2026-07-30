# Tutorial 3: mathematical model

This is the only tutorial containing workflow equations. Other documentation
links here instead of repeating mathematical definitions.

## 1. Condition-aware Pando model

For one broad cell type, target gene \(g\), condition \(c\), candidate edge \(e\)
and paired cell \(i\), the predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

All conditions of the cell type share the same candidate supergraph and
coordinate system.

### Equal-condition transform

For \(K\) conditions,

\[
\mu_e=\frac{1}{K}\sum_c\bar x_{e,c},
\qquad
s_e^2=\frac{1}{K}\sum_c\frac{1}{n_c}
\sum_{i\in c}(x_{e,i}-\bar x_{e,c})^2,
\]

\[
z_{e,i}=\frac{x_{e,i}-\mu_e}{s_e}.
\]

The target response uses the analogous equal-condition centre and
within-condition variance scale. Every outer fold estimates these transforms
from training cells only.

### Sparse-group multitask objective

Let \(B\) contain one coefficient column per condition. For fixed \(\lambda\),
\(\alpha\), and \(ho=condition\_mix\), Pando minimizes

\[
\sum_c\frac{1}{2n_c}\lVert y_c-a_c-X_c\beta_c\rVert_2^2
+\frac{\lambda(1-\alpha)}{2}\lVert B\rVert_F^2
+\lambda\alpha\left[(1-\rho)\sum_e\lVert B_{e,\cdot}\rVert_2
+\rho\sum_{e,c}|\beta_{e,c}|\right].
\]

The row-wise term couples selection across conditions; the entry-wise term
permits condition-specific zeros. A support-constrained common-metric ridge
refit produces the final condition coefficients.

## 2. Estimability and projectable structural zeros

For outer fold \(k\), define the coefficient-estimability indicator

\[
a^{(-k)}_{e,c}=
\mathbf 1\left\{Var_{train(-k,c)}(x_e)>0\right\}.
\]

A non-estimable candidate remains projectable with deterministic zero
contribution:

\[
s^{(-k)}_{e,c}=1-a^{(-k)}_{e,c},
\qquad
\beta^{(-k)}_{e,c}=NA,
\qquad
x^{proj}_{e,i,c}=0\quad\text{when }s^{(-k)}_{e,c}=1.
\]

Thus coefficient availability and projection availability are distinct:

\[
projection\_support_{e,c}=a_{e,c}\lor s_{e,c}.
\]

An exact-zero predictor, such as a peak closed in every input cell, remains in
the shared candidate supergraph. It receives no fitted coefficient and
contributes zero in every condition. The zero status of held-out cells is
determined from the corresponding training fold only.

## 3. Primary condition-full OOF projection

For held-out cell \(i\), target \(g\), condition \(c\), and outer fold \(k\), the
primary regulatory score is

\[
G^{full}_{i,g,c}=
\sum_{e:\,target(e)=g}
a^{(-k)}_{e,c}\,\beta^{(-k)}_{e,c}z^{(-k)}_{i,e}.
\]

Every non-estimable edge side contributes zero. Therefore unilateral edges
contribute only in the estimable condition, and bilaterally non-estimable edges
contribute zero in both conditions.

For a requested comparison set \(C^*\), the common-support component is

\[
m^{(-k)}_e=
\prod_{c\in C^*}a^{(-k)}_{e,c},
\]

\[
G^{common}_{i,g,c}=
\sum_{e:\,target(e)=g}
m^{(-k)}_e\,
\beta^{(-k)}_{e,c}z^{(-k)}_{i,e}.
\]

The condition-unique projection component is

\[
G^{unique}_{i,g,c}=
G^{full}_{i,g,c}-G^{common}_{i,g,c}.
\]

`condition_full_oof` is the primary penalty input. Common support is retained as
a decomposition, not as the primary route.

For metacell \(u\) with exact membership set \(M_u\),

\[
G_{u,g,c}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c}.
\]

Projection is always computed before metacell aggregation.

## 4. Reliability and calibration

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

For target \(g\) and cell type \(t\), RegCompass estimates one robust scale from
the primary condition-full metacell projection:

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G^{full}_{g,t})}{1.349},
MAD_{1.4826}(G^{full}_{g,t}),
\sqrt{mean((G^{full}_{g,t})^2)},
10^{-6}
\right).
\]

The bounded signed regulatory modifier is

\[
R_{g,c,u}=q_g\tanh\left(
\frac{G^{full}_{u,g,c}}{\sigma_{g,t}}
\right).
\]

## 5. RNA and multiome gene support

Metacell RNA expression is converted to bounded support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

The regulatory modifier acts on RNA-support odds:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{R_{g,u}}}.
\]

A non-finite target modifier falls back to the neutral value \(R=0\), which is
exactly RNA-only support. RNA support equal to 0 or 1 remains unchanged.

## 6. GPR aggregation and reaction penalty

For GPR AND group \(A_{r,j}\), the canonical limiting-subunit rule is

\[
Q_{r,j,u}=\min_{g\in A_{r,j}}C^{MO}_{g,u}.
\]

Isozyme OR branches are additive:

\[
E_{r,u}=\sum_jQ_{r,j,u}.
\]

Reaction support is converted to the LP penalty

\[
p_{r,u}=\frac{1}{1+\log_2(1+E_{r,u})}.
\]

## 7. Shared metabolic model

For each medium, RegCompass applies one global FASTCORE completion to the merged
reaction catalogue and reuses the resulting stoichiometric model, bounds,
reaction order and target definitions for every condition and metacell.

For target reaction \(r\) and direction \(d\), the first LP computes

\[
v^{max}_{r,d}=\max_v d v_r
\]

subject to

\[
Sv=0,\qquad l\le v\le u.
\]

The second LP minimizes the network-wide penalty while maintaining a fraction
\(\omega\) of the target maximum:

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
target direction; it is not measured flux.

## 8. Condition statistics

RegCompass transforms normalized penalty to

\[
S_{r,d,u}=-\log(\widetilde P_{r,d,u}+\epsilon).
\]

Pairwise tests use Wilcoxon rank-sum statistics. Analyses with at least three
conditions may use a Kruskal-Wallis omnibus test. Metacells are the statistical
units, so reported P values describe within-dataset separation and are not
donor-level biological-replicate inference.

## 9. Canonical scope

The primary schema contains the condition-full, common-support and RNA-only
routes on one shared structural model. It does not calculate or persist:

- depth matching;
- common-depth restriction;
- alpha sensitivity;
- zero-support sensitivity;
- link-saturation propagation.

Public API: [functions.md](functions.md).
