# RegCompass mathematical specification

This file is the mathematical specification for the RegCompass workflow.

## 1. Condition-comparable regulatory model

For one broad cell type, target gene \(g\), pooled cells, and conditions \(c\), let \(E_g^{pool}\) be the pooled candidate set and \(E_{g,c}\) the candidate set found within condition \(c\). The frozen dictionary is the union of complete TF–peak–target triples:

\[
E_g^{\cup}=E_g^{pool}\cup\bigcup_c E_{g,c}.
\]

For edge \(e\) and paired cell \(i\), the unscaled predictor is

\[
x_{e,i}=RNA_{TF(e),i}\,ATAC_{peak(e),i}.
\]

Every eligible condition uses the same ordered dictionary and fits

\[
y_{g,i}=\alpha_{g,c}+\sum_{e\in E_g^{\cup}}\beta_{e,g,c}x_{e,i}+\varepsilon_{g,i}.
\]

The implementation uses a Gaussian identity GLM with an intercept, `interaction_term = ":"`, and `scale = FALSE`. Pooled coefficients are not used to calibrate condition coefficients.

Let \(a_{e,g,c}\) indicate that an edge coefficient is estimable and let \(m_{g,c}\) indicate that the complete target-level regression has `fit_status == "ok"`. The adjusted-P threshold is

\[
q=0.05.
\]

An edge is active for RegCompass regulatory penalty projection only when the target regression is full-rank under the fitted dictionary, the coefficient itself is estimable, the estimate is finite, and the BH-adjusted P value is below 0.05:

\[
s_{e,g,c}=\mathbf{1}\left\{
 m_{g,c}=1
 \ \land\ a_{e,g,c}=1
 \ \land\ padj_{e,g,c}<q
 \ \land\ \widehat\beta_{e,g,c}\text{ is finite}
\right\}.
\]

No additional post-fit absolute-correlation or absolute-effect-size threshold is applied. Candidate TF and peak correlation filters act upstream during Pando candidate discovery, not as a second coefficient gate.

The effect used for downstream projection is

\[
\theta_{e,g,c}=
\begin{cases}
\widehat\beta_{e,g,c}, & s_{e,g,c}=1,\\
0, & s_{e,g,c}=0.
\end{cases}
\]

A `rank_deficient` target remains in the complete coefficient and fit-diagnostic tables for audit, but every edge belonging to that target has zero realized RegCompass penalty contribution even when an individual coefficient is finite and its P value would otherwise pass BH. This policy avoids attributing an edge-specific regulatory effect when the complete frozen-dictionary coefficient vector is not uniquely identifiable in that condition. Other non-`ok` target statuses, including insufficient residual degrees of freedom and failed/non-finite fits, are likewise excluded from penalty projection.

The Pando source object may retain its original GLM significance fields. RegCompass records the target fit status on the gated coefficient table and applies the stricter `fit_status == "ok"` rule before paired-cell regulatory projection and before active-edge assembly.

## 2. Paired-cell projection and metacell aggregation

For cell \(i\) in condition \(c\), the target-level regulatory score is

\[
G_{i,g,c}=\sum_{e\in E_g^{\cup}}\theta_{e,g,c}x_{e,i}.
\]

For metacell \(u\) with membership set \(M_u\),

\[
G_{u,g}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c(i)}.
\]

For target \(g\) and cell type \(t\), the calibration scale is

\[
\sigma_{g,t}=\max\left(\frac{IQR(G_{g,t})}{1.349},MAD_{1.4826}(G_{g,t}),\sqrt{mean(G_{g,t}^{2})},10^{-6}\right).
\]

The bounded regulatory modifier is

\[
R_{g,u}=q_{g,u}\tanh\left(\frac{G_{g,u}}{\sigma_{g,t(u)}}\right).
\]

Unavailable regulatory evidence uses the neutral RNA-only route downstream.

## 3. RNA and multiome gene support

Metacell RNA support is

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Regulatory evidence modifies the support odds:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{R_{g,u}}}{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{R_{g,u}}}.
\]

## 4. GPR aggregation and COMPASS reaction cost

For an AND branch \(A_{r,j}\),

\[
Q_{r,j,u}=\min_{g\in A_{r,j}}C^{MO}_{g,u}.
\]

Isozyme branches are additive:

\[
E_{r,u}=\sum_jQ_{r,j,u}.
\]

Reaction support is converted to the COMPASS-style LP cost:

\[
p_{r,u}=\frac{1}{1+\log_2(1+\max(E_{r,u},0))}.
\]

If reaction expression is unavailable, RegCompass sets \(E_{r,u}=0\), so

\[
p_{r,u}=1.
\]

The structural roles `exchange`, `demand`, `sink`, and `artificial_support` also receive the maximum COMPASS cost

\[
p_{r,u}=1,
\]

rather than a separate fixed cost outside the COMPASS range.

## 5. Default CORDA2 structural model

For cell type \(t\) and condition-specific biological reaction sets \(B_{t,c}\), the structural catalogue is

\[
B_t=\bigcup_c B_{t,c}.
\]

Layer 1 reaction evidence is summarized within condition and unioned within cell type, then mapped to the CORDA2 confidence groups

\[
HC_{t},\ MC_{t},\ NC_{t},\ OT_{t}.
\]

For medium \(m\), exchange bounds are first applied to the complete parent GEM. The default cell-type structural model is

\[
\mathcal{G}_{t,m}=
CORDA2(HC_t,MC_t,NC_t,OT_t;S,l_m,u_m).
\]

The parent is passed directly to the original MATLAB CORDA2 state machine without FASTCC pre-pruning. Reversible reactions are represented by non-negative directional variables during reconstruction. When selected directions are merged back into reactions, excluded directions are fixed to zero and retained directions recover the corresponding parent bounds. In particular, a retained irreversible reaction with positive parent lower bound \(l_{m,r}>0\) retains

\[
l^{final}_{t,m,r}=l_{m,r}>0.
\]

FASTCORE is an explicit supplementary reconstruction selected with `model_completion = "fastcore"`. The full-GEM route is a supplementary mode that keeps the complete medium-constrained parent and performs no context-specific reconstruction.

## 6. Directional COMPASS-like LP

For target reaction \(r\) and direction \(d\),

\[
v^{max}_{r,d,t,m}=\max_v d\,v_r
\]

subject to

\[
S_{t,m}v=0,\qquad l_{t,m}\le v\le u_{t,m}.
\]

For metacell \(u\), RegCompass solves

\[
P^{*}_{r,d,u,m}=\min_{v,z}\sum_jp_{j,u}z_j
\]

subject to

\[
S_{t,m}v=0,\qquad -z_j\le v_j\le z_j,\qquad d\,v_r\ge\omega v^{max}_{r,d,t,m}.
\]

The normalized penalty is

\[
\widetilde P_{r,d,u,m}=\frac{P^{*}_{r,d,u,m}}{\omega v^{max}_{r,d,t,m}}.
\]

Lower normalized penalty indicates stronger model-constrained support for the target direction; it is not a measured flux.

## 7. Condition statistics

The comparison score is

\[
S_{r,d,u,m}=-\log(\widetilde P_{r,d,u,m}+\epsilon).
\]

Condition comparisons are performed within one cell type, target reaction, direction, and medium. Wilcoxon rank-sum tests are used for pairwise comparisons; Kruskal–Wallis tests may be used for analyses with at least three conditions. Metacells are within-dataset statistical units, not donor-level biological replicates.

## 8. Inference scope

Condition-GRN P values are conditional on the frozen candidate dictionary and do not include selective-inference correction for candidate discovery. Standard-Pando cell types do not receive manufactured condition coefficients.
