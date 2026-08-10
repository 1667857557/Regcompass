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

Pando coefficient magnitudes are retained in the GRN fit for inference diagnostics and audit, but **absolute coefficient magnitude is not a RegCompass penalty weight**. The downstream signed edge direction is

\[
d_{e,g,c}=
\begin{cases}
\operatorname{sign}(\widehat\beta_{e,g,c}), & s_{e,g,c}=1,\\
0, & s_{e,g,c}=0.
\end{cases}
\]

Thus two active positive edges with estimates 0.2 and 2.0 have the same RegCompass directional weight \(+1\); two active negative edges with estimates -0.2 and -2.0 have the same directional weight \(-1\). Estimate magnitude cannot change the Layer-1 penalty projection after the significance/estimability gate.

A `rank_deficient` target remains in the complete coefficient and fit-diagnostic tables for audit, but every edge belonging to that target has zero realized RegCompass penalty contribution even when an individual coefficient is finite and its P value would otherwise pass BH. This policy avoids attributing an edge-specific regulatory direction when the complete frozen-dictionary coefficient vector is not uniquely identifiable in that condition. Other non-`ok` target statuses, including insufficient residual degrees of freedom and failed/non-finite fits, are likewise excluded from penalty projection.

The Pando source object may retain its original GLM coefficient and significance fields. RegCompass records the stricter target fit status and applies `fit_status == "ok"` before signed paired-cell regulatory projection and before active-edge assembly.

## 2. Sign-only paired-cell projection and metacell aggregation

Removing \(|\widehat\beta|\) also removes the coefficient's automatic compensation for arbitrary predictor units. RegCompass therefore calibrates each TF-by-ATAC predictor using only its own distribution over **all fitted cells of the same broad cell type**, pooled across conditions. This preserves condition comparability while preventing one edge from dominating only because its normalized RNA or ATAC feature has a larger numerical scale.

For edge \(e\) in broad cell type \(t\), define

\[
\lambda_{e,t}=\max\left(
\frac{IQR(x_{e,t})}{1.349},
MAD_{1.4826}(x_{e,t}),
\sqrt{mean(x_{e,t}^{2})},
10^{-6}
\right),
\]

with a neutral scale of 1 when the predictor is identically zero. The bounded paired-cell edge activity is

\[
a_{e,i}=\tanh\left(\frac{\max(x_{e,i},0)}{\lambda_{e,t(i)}}\right).
\]

Because both \(x\) and \(\lambda\) multiply by the same positive constant under a change of feature units, \(a_{e,i}\) is invariant to positive rescaling of the raw TF-by-ATAC predictor. It also retains paired-cell co-occurrence because TF and ATAC are multiplied before SuperCell aggregation.

Let the number of active edges for target \(g\) in condition \(c\) be

\[
N_{g,c}=\sum_{e\in E_g^{\cup}}\mathbf{1}\{d_{e,g,c}\neq0\}.
\]

For \(N_{g,c}>0\), the target-level sign-only regulatory score for cell \(i\) is the **mean signed active-edge activity**

\[
G_{i,g,c}=\frac{1}{N_{g,c}}
\sum_{e\in E_g^{\cup}}d_{e,g,c}\,a_{e,i}.
\]

This normalization is essential after coefficient magnitudes are removed: otherwise a target or condition with more significant Pando edges would receive a larger regulatory score solely because of network degree. Because \(|d_{e,g,c}a_{e,i}|\le1\), the degree-normalized cell score satisfies

\[
-1\le G_{i,g,c}\le1.
\]

For metacell \(u\) with exact membership set \(M_u\),

\[
G_{u,g}=\frac{1}{|M_u|}\sum_{i\in M_u}G_{i,g,c(i)}.
\]

Stage 2 remains condition-pure, so each final metacell receives the edge directions and active-edge degree fitted for its own condition. The membership operation is unchanged from the previous Layer-1 contract; only the edge magnitude definition changes from coefficient-weighted to degree-controlled sign-only, self-scaled paired-cell activity. If a target has no active edge in a metacell's condition, regulatory evidence is unavailable and the downstream RNA-only fallback remains neutral.

For target \(g\) and cell type \(t\), the existing target-level calibration scale remains

\[
\sigma_{g,t}=\max\left(\frac{IQR(G_{g,t})}{1.349},MAD_{1.4826}(G_{g,t}),\sqrt{mean(G_{g,t}^{2})},10^{-6}\right).
\]

The bounded regulatory modifier remains

\[
R_{g,u}=q_{g,u}\tanh\left(\frac{G_{g,u}}{\sigma_{g,t(u)}}\right).
\]

The reliability term \(q_{g,u}\) retains the existing routing semantics: condition-Pando uses reliability 1 for targets with at least one mapped active edge in that condition, whereas standard Pando retains its target-level fit reliability. Unavailable regulatory evidence uses the neutral RNA-only route downstream.

For the one-effective-condition standard-Pando route, the same sign-only paired-cell activity and active-edge-degree normalization are used after the standard Pando adjusted-P/estimability gate. Standard-Pando coefficient magnitudes therefore also remain audit information rather than penalty magnitudes.

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
