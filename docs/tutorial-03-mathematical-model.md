# Tutorial 3: mathematical model

This tutorial collects the mathematical definitions used by the complete
RegCompass workflow. The condition-GRN sections reflect the common-dictionary
Pando model; the RNA-support, GPR, metabolic LP, and condition-statistics
sections are retained and indexed by the structural cell-type boundary.

## 1. Common-dictionary condition GRN

Consider one broad cell type \(t\), target gene \(g\), condition \(c\), candidate
edge \(e\), and paired cell \(i\). An edge is an exact TF–peak–target triple

\[
e=(TF(e),peak(e),g).
\]

Candidate discovery is run on all eligible cells of cell type \(t\) and
separately within each condition. Let \(D_{g,t}^{pool}\) denote the pooled
candidate set and \(D_{g,t,c}\) the candidate set from condition \(c\). The
common reference is the exact triple union

\[
D_{g,t}^{\cup}=D_{g,t}^{pool}\cup\bigcup_c D_{g,t,c}.
\]

The union is performed on complete triples. TFs, peaks, and targets are not
recombined by Cartesian product. The resulting cell-type- and target-specific
dictionary is frozen before coefficient estimation.

For every edge in the frozen dictionary, the unscaled interaction predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

RNA and ATAC preprocessing is completed before condition splitting, so every
condition of the same cell type uses the same feature definitions and numerical
units. For each condition, Pando fits

\[
y_{g,i}=a_{g,t,c}+\sum_{e\in D_{g,t}^{\cup}}\beta_{e,t,c}x_{e,i}
+\varepsilon_{g,i},
\qquad i\in(t,c),
\]

with a Gaussian identity GLM, an intercept, `interaction_term = ":"`, and
`scale = FALSE`. The condition effect is the fitted coefficient
\(\widehat\beta_{e,t,c}\). A pooled coefficient is not used to centre, scale, or
calibrate the condition coefficient.

The candidate-discovery correlations (`tf_cor` and `peak_cor`) are not rerun
after the dictionary is frozen. Consequently, all conditions of cell type \(t\)
are fitted against the same ordered predictor dictionary even when significant
edge sets or coefficient directions differ.

## 2. Estimability, uncertainty, and penalty-edge selection

An edge is estimable in condition \(c\) only when its predictor has usable
variance and its coefficient is identifiable in the condition-specific design
matrix. Define

\[
a_{e,t,c}=\mathbf 1\{\widehat\beta_{e,t,c}\text{ is estimable}\}.
\]

A zero-variance, aliased, non-finite, or insufficient-residual-df edge remains

\[
\widehat\beta_{e,t,c}=NA,\qquad a_{e,t,c}=0.
\]

It is not converted into a fitted zero and is not interpreted as evidence for
absence of regulation.

For estimable coefficients, Pando retains the standard error, test statistic,
raw P value, and adjusted P value. Benjamini-Hochberg adjustment is performed
within the fitted condition network. The penalty-selection indicator is

\[
s_{e,t,c}=\mathbf 1\left\{
 a_{e,t,c}=1\ \land\ padj_{e,t,c}<0.05
\right\}.
\]

The coefficient passed to RegCompass is

\[
\theta_{e,t,c}=
\begin{cases}
\widehat\beta_{e,t,c}, & s_{e,t,c}=1,\\
0, & s_{e,t,c}=0.
\end{cases}
\]

The complete coefficient table is retained for audit. A non-significant edge
therefore remains distinguishable from an unavailable edge and from an
estimable coefficient that is numerically close to zero.

These ordinary GLM P values are conditional on the selected frozen dictionary;
they do not include a selective-inference correction for the candidate-discovery
step.

## 3. Cell-first regulatory projection and metacell aggregation

For a cell of type \(t\) in condition \(c\), the regulatory contribution to
target \(g\) is

\[
G_{i,g,t,c}=\sum_{e\in D_{g,t}^{\cup}}\theta_{e,t,c}x_{e,i}
=\sum_{e\in D_{g,t}^{\cup}}
\theta_{e,t,c}RNA_{TF(e),i}ATAC_{peak(e),i}.
\]

The coefficient is always selected from the cell's own condition. RegCompass
computes this quantity on paired cells before aggregation; it does not replace
\(RNA\times ATAC\) with a product of metacell averages.

For metacell \(u\) with exact SuperCell membership set \(M_u\), cell type
\(t(u)\), and condition \(c(u)\),

\[
G_{u,g}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,t(u),c(u)}.
\]

A condition–target combination without a significant estimable edge has no
finite regulatory modifier and uses the neutral RNA-only fallback downstream.
It is recorded explicitly in projection coverage rather than represented by a
manufactured coefficient.

Several output fields retain historical names containing `_oof`, `common`, or
`condition_unique` for API compatibility. In the current model:

\[
G^{primary}=G^{common}=G,
\qquad
G^{condition\_unique}=0.
\]

These aliases do not imply outer-fold fitting, a shared-slope model, or a
common-support decomposition.

## 4. Projection calibration and regulatory modifier

For target \(g\) and broad cell type \(t\), RegCompass estimates a robust scale
from finite primary metacell projections:

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},
10^{-6}
\right).
\]

For a condition–target combination with at least one significant estimable
edge, reliability is one. If no such edge exists, reliability is unavailable:

\[
q_{g,t,c}=\begin{cases}
1, & \exists e:s_{e,t,c}=1,\\
NA, & \text{otherwise}.
\end{cases}
\]

The bounded signed regulatory modifier is

\[
R_{g,u}=q_{g,t(u),c(u)}\tanh\left(
\frac{G_{u,g}}{\sigma_{g,t(u)}}
\right).
\]

A non-finite modifier is handled as unavailable regulatory evidence and is
replaced by the neutral value \(R=0\) only at the RNA/multiome integration step.

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

## 7. Cell-type-specific union GEM and FASTCORE

Let \(B_{t,c}\) be the biological reaction set produced for cell type \(t\) and
condition \(c\). Stage 3 unions conditions only within cell type:

\[
B_t=\bigcup_c B_{t,c}.
\]

There is no cross-cell-type union \(\bigcup_t B_t\). Different cell types may
therefore retain different core reactions and different biological expansion
sets.

For medium \(m\), RegCompass applies the medium bounds to the parent GEM and
runs FASTCORE independently for each cell type:

\[
\mathcal G_{t,m}=FASTCORE(B_t;S,l_m,u_m).
\]

Thus, the structural cache key is \((t,m)\), not only \(m\). The resulting model
file, checksum, reaction order, support reactions and bounds are specific to
that cell type and medium. Conditions of cell type \(t\) share
\(\mathcal G_{t,m}\); a metacell from another cell type cannot use it.

For target reaction \(r\), direction \(d\), cell type \(t\), and medium \(m\),
the first LP computes

\[
v^{max}_{r,d,t,m}=\max_v d v_r
\]

subject to the matching structural model

\[
S_{t,m}v=0,\qquad l_{t,m}\le v\le u_{t,m}.
\]

The second LP for metacell \(u\), where \(t(u)=t\), minimizes the network-wide
penalty while maintaining a fraction \(\omega\) of the target maximum:

\[
P^*_{r,d,u,m}=\min_{v,z}\sum_jp_{j,u}z_j
\]

subject to

\[
S_{t,m}v=0,\qquad -z_j\le v_j\le z_j,
\qquad d v_r\ge\omega v^{max}_{r,d,t,m}.
\]

The normalized penalty is

\[
\widetilde P_{r,d,u,m}=\frac{P^*_{r,d,u,m}}
{\omega v^{max}_{r,d,t,m}}.
\]

Directional \(v^{max}\) is solved once per cell-type model and target direction,
then reused across conditions and metacells of the same cell type. Lower
normalized penalty indicates stronger network-constrained support for the target
direction; it is not measured flux.

## 8. Condition statistics

RegCompass transforms normalized penalty to

\[
S_{r,d,u,m}=-\log(\widetilde P_{r,d,u,m}+\epsilon).
\]

A condition contrast fixes cell type \(t\), reaction \(r\), direction \(d\), and
medium \(m\). Both conditions therefore use the same
\(\mathcal G_{t,m}\), bounds and \(v^{max}_{r,d,t,m}\). Rows belonging to a
different cell type are excluded; they are not treated as missing observations
on a shared global model.

Pairwise tests use Wilcoxon rank-sum statistics. Analyses with at least three
conditions may use a Kruskal-Wallis omnibus test. Metacells are the statistical
units, so reported P values describe within-dataset separation and are not
donor-level biological-replicate inference.

Condition-GRN coefficients may also be compared descriptively on the common
edge dictionary within one cell type:

\[
\Delta\beta_{e,t,c_1,c_2}=
\widehat\beta_{e,t,c_1}-\widehat\beta_{e,t,c_2}.
\]

Such comparisons require checking estimability and uncertainty in both
conditions. A non-significant or unavailable coefficient must not be interpreted
as a biological zero.

## 9. Canonical scope

The primary schema contains the BH-filtered fixed-dictionary condition route and
an RNA-only fallback/control. Structural metabolic models are shared across
conditions only within one cell type and medium. Historical `common` fields are
aliases of the primary route and historical `condition_unique` fields are zero
compatibility matrices.

The workflow does not calculate or persist:

- nested or outer-fold condition-GRN estimates;
- sparse-group condition paths;
- structural-zero condition coefficients;
- pooled-coefficient calibration;
- cross-cell-type meta-module unions;
- cross-cell-type union GEMs;
- global FASTCORE completion across cell types;
- depth matching;
- common-depth restriction;
- alpha sensitivity;
- zero-support sensitivity;
- link-saturation propagation.

Public API: [functions.md](functions.md).
