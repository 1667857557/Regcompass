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

RegCompass constrains the realized modifier to

\[
-1\le R_{g,u}\le1.
\]

Unavailable regulatory evidence uses the neutral RNA-only route downstream, equivalent to \(R_{g,u}=0\).

## 3. Two distinct Layer 1 evidence scales

Layer 1 intentionally exposes two different numerical representations. They must not be interchanged because they answer different questions.

### 3.1 Quantitative COMPASS penalty path

Let \(\widehat{CPM}_{g,u}\) denote the empirical-Bayes latent metacell CPM produced by the RNA model. The quantitative RNA expression used for the LP penalty is

\[
X^{RNA}_{g,u}=\widehat{CPM}_{g,u}.
\]

Pando modifies this unbounded non-negative expression multiplicatively:

\[
\boxed{
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}}
}
\]

so that

\[
\frac12X^{RNA}_{g,u}\le X^{MO}_{g,u}\le2X^{RNA}_{g,u}.
\]

This is the quantitative regulatory correction. It preserves the latent-CPM dynamic range instead of first mapping every gene into \([0,1]\). It also preserves the zero condition:

\[
X^{RNA}_{g,u}=0\Rightarrow X^{MO}_{g,u}=0.
\]

For an AND branch \(A_{r,j}\), the canonical RegCompass workflow uses the configured GPR AND operator, `min` by default:

\[
Q^{quant}_{r,j,u}=\min_{g\in A_{r,j}}X^{MO}_{g,u}.
\]

Isozyme branches remain additive, matching the existing COMPASS-style OR policy:

\[
E^{quant}_{r,u}=\sum_jQ^{quant}_{r,j,u}.
\]

The reaction value entering the LP cost is therefore unbounded above when expression is high. RegCompass then applies the COMPASS-shaped cost exactly once at reaction level:

\[
\boxed{
p_{r,u}=\frac{1}{1+\log_2(1+\max(E^{quant}_{r,u},0))}
}
\]

If quantitative reaction expression is unavailable, RegCompass sets \(E^{quant}_{r,u}=0\), hence

\[
p_{r,u}=1.
\]

The structural roles `exchange`, `demand`, `sink`, and `artificial_support` also receive

\[
p_{r,u}=1.
\]

This restores a COMPASS-like expression dynamic range: for a single-gene, single-isozyme reaction, increasing expression can drive \(p_{r,u}\) toward zero rather than imposing a floor near 0.5. RegCompass still differs from the original COMPASS implementation in other design choices, including the default GPR AND rule (`min` rather than the original default `mean`) and the empirical-Bayes metacell expression model.

### 3.2 Bounded structural-confidence path

A separate bounded representation is retained for CORDA2 and other structural-confidence decisions. Let

\[
L_{g,u}=\ln(1+\widehat{CPM}_{g,u}).
\]

The bounded RNA support is

\[
C^{RNA}_{g,u}=\frac{L_{g,u}}{L_{g,u}+h}.
\]

Regulatory evidence modifies the support odds:

\[
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{R_{g,u}}}{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{R_{g,u}}}.
\]

Thus

\[
0\le C^{RNA}_{g,u},C^{MO}_{g,u}\le1.
\]

Bounded reaction structural support is then computed with the same GPR topology:

\[
Q^{struct}_{r,j,u}=\min_{g\in A_{r,j}}C^{MO}_{g,u},
\qquad
E^{struct}_{r,u}=\sum_jQ^{struct}_{r,j,u}.
\]

The bounded representation is **not** passed to the quantitative LP penalty in Layer 1 schema v5. It is retained to stabilize structural evidence and preserve the previous CORDA2 confidence scale.

For backward compatibility, Layer 1 continues to expose the bounded matrices as `reaction_expression` and `reaction_expression_rna_only`, with explicit aliases `reaction_structural_support` and `reaction_structural_support_rna_only`. The quantitative LP matrices are `reaction_expression_quantitative` and `reaction_expression_quantitative_rna_only`. Layer 2 explicitly prefers the quantitative matrices when schema v5 is present.

## 4. Why the two paths are separated

The previous formulation used bounded gene support for both structural confidence and quantitative LP cost:

\[
\widehat{CPM}\rightarrow\ln(1+\widehat{CPM})
\rightarrow\frac{L}{L+h}\rightarrow GPR
\rightarrow\frac{1}{1+\log_2(1+E)}.
\]

For a single-gene, single-isozyme reaction this forced \(E<1\), so even arbitrarily high RNA expression could not reduce the expression penalty below approximately 0.5. It also attenuated the realized Pando effect at high expression because the regulatory correction acted inside a saturating \([0,1]\) support scale.

Layer 1 schema v5 therefore separates the estimands:

\[
\boxed{
\widehat{CPM}\times2^R\rightarrow GPR\rightarrow COMPASS\ penalty
}
\]

for quantitative metabolic scoring, and

\[
\boxed{
\frac{\ln(1+\widehat{CPM})}{\ln(1+\widehat{CPM})+h}
\xrightarrow{\text{bounded odds }2^R}
structural\ support
}
\]

for CORDA2/structural confidence.

This change does not alter Pando fitting, metacell construction, the frozen edge dictionary, rank-deficiency gating, CORDA2 reconstruction semantics, directional \(V_{max}\), or the LP constraints. It changes only which Layer 1 numerical evidence scale is supplied to the LP objective.

## 5. Default CORDA2 structural model

For cell type \(t\) and condition-specific biological reaction sets \(B_{t,c}\), the structural catalogue is

\[
B_t=\bigcup_c B_{t,c}.
\]

CORDA2 uses the bounded structural reaction evidence \(E^{struct}\), not \(E^{quant}\). Layer 1 structural evidence is summarized within condition and unioned within cell type, then mapped to the CORDA2 confidence groups

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

## 7. RNA-only interpretation control

The RNA-only comparison reuses the exact same structural model, media, bounds, target directions, and directional \(V_{max}\) as the multiome route. Only the quantitative objective coefficients differ:

\[
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}}
\]

versus

\[
X^{RNA}_{g,u}.
\]

Therefore `regulatory_penalty_delta` isolates the effect of the Pando regulatory modifier on the quantitative COMPASS objective **conditional on the already constructed shared structural model**. It is not an end-to-end RNA-only reconstruction control.

## 8. Condition statistics

The comparison score is

\[
S_{r,d,u,m}=-\log(\widetilde P_{r,d,u,m}+\epsilon).
\]

Condition comparisons are performed within one cell type, target reaction, direction, and medium. Wilcoxon rank-sum tests are used for pairwise comparisons; Kruskal–Wallis tests may be used for analyses with at least three conditions. Metacells are within-dataset statistical units, not donor-level biological replicates.

## 9. Inference scope

Condition-GRN P values are conditional on the frozen candidate dictionary and do not include selective-inference correction for candidate discovery. Standard-Pando cell types do not receive manufactured condition coefficients.

Layer 1 schema v5 is intentionally incompatible with old v4 artifacts for downstream validation. Existing `step_layer1.rds` objects must be regenerated before Layer 2 so that the quantitative reaction-expression matrices are present and the LP cannot silently fall back to the previous bounded penalty scale.
