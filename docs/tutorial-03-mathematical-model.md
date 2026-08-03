# Tutorial 3: mathematical model

This tutorial collects the mathematical definitions used by the complete
RegCompass workflow. The condition-GRN sections reflect the common-dictionary
Pando model; the RNA-support, GPR, metabolic LP, and condition-statistics
sections are retained because those parts of the workflow are unchanged.

## 1. Common-dictionary condition GRN

Consider one broad cell type, target gene \(g\), condition \(c\), candidate edge
\(e\), and paired cell \(i\). An edge is an exact TF–peak–target triple

\[
e=(TF(e),peak(e),g).
\]

Candidate discovery is run on all eligible cells of the cell type and separately
within each condition. Let \(D_g^{global}\) denote the pooled candidate set and
\(D_{g,c}\) the candidate set from condition \(c\). The common reference is the
exact triple union

\[
D_g^{\cup}=D_g^{global}\cup\bigcup_c D_{g,c}.
\]

The union is performed on complete triples. TFs, peaks, and targets are not
recombined by Cartesian product. The resulting target-specific dictionary is
frozen before coefficient estimation.

For every edge in the frozen dictionary, the unscaled interaction predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

RNA and ATAC preprocessing is completed before condition splitting, so every
condition uses the same feature definitions and numerical units. For each
condition, Pando fits

\[
y_{g,i}=a_{g,c}+\sum_{e\in D_g^{\cup}}\beta_{e,c}x_{e,i}
+\varepsilon_{g,i},
\qquad i\in c,
\]

with a Gaussian identity GLM, an intercept, `interaction_term = ":"`, and
`scale = FALSE`. The condition effect is the fitted coefficient
\(\widehat\beta_{e,c}\). A pooled coefficient is not used to centre, scale, or
calibrate the condition coefficient.

The candidate-discovery correlations (`tf_cor` and `peak_cor`) are not rerun
after the dictionary is frozen. Consequently, all conditions are fitted against
the same ordered predictor dictionary even when the significant edge sets or
coefficient directions differ.

## 2. Estimability, uncertainty, and penalty-edge selection

An edge is estimable in condition \(c\) only when its predictor has usable
variance and its coefficient is identifiable in the condition-specific design
matrix. Define

\[
a_{e,c}=\mathbf 1\{\widehat\beta_{e,c}\text{ is estimable}\}.
\]

A zero-variance, aliased, non-finite, or insufficient-residual-df edge remains

\[
\widehat\beta_{e,c}=NA,\qquad a_{e,c}=0.
\]

It is not converted into a fitted zero and is not interpreted as evidence for
absence of regulation.

For estimable coefficients, Pando retains the standard error, test statistic,
raw P value, and adjusted P value. Benjamini-Hochberg adjustment is performed
within the fitted condition network. The penalty-selection indicator is

\[
s_{e,c}=\mathbf 1\left\{
 a_{e,c}=1\ \land\ padj_{e,c}<0.05
\right\}.
\]

The coefficient passed to RegCompass is

\[
\theta_{e,c}=
\begin{cases}
\widehat\beta_{e,c}, & s_{e,c}=1,\\
0, & s_{e,c}=0.
\end{cases}
\]

The complete coefficient table is retained for audit. A non-significant edge
therefore remains distinguishable from an unavailable edge and from an
estimable coefficient that is numerically close to zero.

These ordinary GLM P values are conditional on the selected frozen dictionary;
they do not include a selective-inference correction for the candidate-discovery
step.

## 3. Cell-first regulatory projection and metacell aggregation

For a cell in condition \(c\), the regulatory contribution to target \(g\) is

\[
G_{i,g,c}=\sum_{e\in D_g^{\cup}}\theta_{e,c}x_{e,i}
=\sum_{e\in D_g^{\cup}}
\theta_{e,c}RNA_{TF(e),i}ATAC_{peak(e),i}.
\]

The coefficient is always selected from the cell's own condition. RegCompass
computes this quantity on paired cells before aggregation; it does not replace
\(RNA\times ATAC\) with a product of metacell averages.

For metacell \(u\) with exact SuperCell membership set \(M_u\),

\[
G_{u,g,c}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c}.
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
from the finite primary metacell projections:

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
q_{g,c}=\begin{cases}
1, & \exists e:s_{e,c}=1,\\
NA, & \text{otherwise}.
\end{cases}
\]

The bounded signed regulatory modifier is

\[
R_{g,c,u}=q_{g,c}\tanh\left(\frac{G_{u,g,c}}{\sigma_{g,t}}\right).
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

## 7. Shared metabolic model

For each medium, RegCompass applies one global FASTCORE completion to the merged
reaction catalogue and reuses the resulting stoichiometric model, bounds,
reaction order, and target definitions for every condition and metacell.

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

Condition-GRN coefficients may also be compared descriptively on the common
edge dictionary:

\[
\Delta\beta_{e,c_1,c_2}=
\widehat\beta_{e,c_1}-\widehat\beta_{e,c_2}.
\]

Such comparisons require checking estimability and uncertainty in both
conditions. A non-significant or unavailable coefficient must not be interpreted
as a biological zero.

## 9. Canonical scope

The primary schema contains the BH-filtered fixed-dictionary condition route and
an RNA-only fallback/control on one shared structural metabolic model. Historical
`common` fields are aliases of the primary route and historical
`condition_unique` fields are zero compatibility matrices.

The workflow does not calculate or persist:

- nested or outer-fold condition-GRN estimates;
- sparse-group condition paths;
- structural-zero condition coefficients;
- pooled-coefficient calibration;
- depth matching;
- common-depth restriction;
- alpha sensitivity;
- zero-support sensitivity;
- link-saturation propagation.

Public API: [functions.md](functions.md).


## Cell-type-specific union GEM

Let `c` denote cell type, `d` condition, and `m` medium. Condition-specific
biological reaction sets are unioned only within cell type,
`B_c = union_d B_{c,d}`. RegCompass then runs FASTCORE independently for every
`(c,m)` pair and obtains `G_{c,m}`. There is no operation of the form
`union_c G_{c,m}`. Consequently, conditions are structurally comparable within
a cell type, while biologically distinct cell types are not forced onto an
artificial cross-cell-type reaction universe.
