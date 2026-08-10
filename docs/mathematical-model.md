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

### 3.1 Quantitative COMPASS penalty path: SuperCell representative RNA state

Let \(Y_{g,i}\) be the raw RNA count for gene \(g\) in original single cell \(i\), and let

\[
L_i=\sum_hY_{h,i}
\]

be that cell's complete RNA library size. RegCompass first computes a **linear**, non-logarithmic per-cell CPM:

\[
x_{g,i}=10^6\frac{Y_{g,i}}{L_i}.
\]

For final SuperCell metacell \(u\), quantitative RNA is the equal-weight mean of its member-cell normalized expression:

\[
\boxed{
X^{RNA}_{g,u}=\frac{1}{|M_u|}\sum_{i\in M_u}x_{g,i}
}
\]

This implements the SuperCell coarse-graining interpretation that a metacell represents the average state of its member cells. It is intentionally different from normalizing the summed metacell counts:

\[
10^6\frac{\sum_{i\in M_u}Y_{g,i}}{\sum_{i\in M_u}L_i}
=
\sum_{i\in M_u}
\frac{L_i}{\sum_jL_j}x_{g,i},
\]

which is a library-size-weighted mean. The quantitative LP path therefore does **not** let a higher-depth cell contribute more merely because more UMIs were captured, and final metacell size does not directly scale enzyme abundance.

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

The Pando modifier is applied after the equal-weight SuperCell RNA aggregation; Pando's regulatory projection itself remains cell-first and is averaged over the same exact membership before the bounded modifier is constructed.

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

The quantitative estimator contains no empirical-Bayes shrinkage and does not borrow RNA abundance across conditions. A zero member-cell count remains zero before averaging; no cell-type prior creates quantitative expression for an unobserved gene. RegCompass still differs from original COMPASS in other design choices, including the default GPR AND rule (`min` rather than the original default `mean`) and the use of SuperCell metacells as statistical units.

### 3.2 Bounded structural-confidence path

A separate bounded representation is retained for CORDA2 and other structural-confidence decisions. In schema v6 this path deliberately retains the pre-existing empirical-Bayes latent metacell CPM; that latent estimator is **structural-only** and is no longer an LP expression input.

Let \(\widehat{CPM}^{struct}_{g,u}\) denote that latent structural-support estimate and

\[
L_{g,u}=\ln(1+\widehat{CPM}^{struct}_{g,u}).
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

The bounded representation is **not** passed to the quantitative LP penalty. It is retained to preserve the established CORDA2 confidence scale while the quantitative RNA estimator is corrected independently.

For backward compatibility, Layer 1 continues to expose the bounded matrices as `reaction_expression` and `reaction_expression_rna_only`, with explicit aliases `reaction_structural_support` and `reaction_structural_support_rna_only`. The quantitative LP matrices are `reaction_expression_quantitative` and `reaction_expression_quantitative_rna_only`. Layer 1 schema v6 additionally saves `rna_metacell_mean_single_cell_cpm`; Layer 2 explicitly prefers the quantitative reaction matrices.

## 4. Why the quantitative RNA estimator changed in schema v6

Schema v5 used

\[
CPM(\text{summed metacell counts})
\rightarrow
\text{cell-type empirical-Bayes shrinkage}
\rightarrow
\widehat{CPM}
\rightarrow
2^R
\rightarrow GPR\rightarrow penalty.
\]

That construction had two undesirable properties for the LP objective. First, `CPM(sum counts)` weights member cells in proportion to their RNA library sizes rather than representing the equal-weight SuperCell state. Second, the empirical-Bayes update shrank metacells toward a shared cell-type mean and could create positive quantitative expression from a pooled zero, thereby attenuating condition differences before metabolic scoring.

Schema v6 separates the estimands more strictly:

\[
\boxed{
Y_{g,i}
\rightarrow
10^6Y_{g,i}/L_i
\xrightarrow{\text{equal SuperCell mean}}
X^{RNA}_{g,u}
\times2^R
\rightarrow GPR
\rightarrow COMPASS\ penalty
}
\]

for quantitative metabolic scoring, and

\[
\boxed{
\widehat{CPM}^{struct}
\rightarrow
\frac{\ln(1+\widehat{CPM}^{struct})}
{\ln(1+\widehat{CPM}^{struct})+h}
\xrightarrow{\text{bounded odds }2^R}
structural\ support
}
\]

for the current CORDA2/structural-confidence route.

This change does not alter Pando fitting, the common edge dictionary, rank-deficiency gating, SuperCell membership construction, CORDA2 reconstruction semantics, directional \(V_{max}\), or the LP constraints. It changes the RNA abundance estimator supplied to the quantitative LP objective and records that estimator explicitly in the Layer 1 contract.

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

## 9. Inference scope and artifact compatibility

Condition-GRN P values are conditional on the frozen candidate dictionary and do not include selective-inference correction for candidate discovery. Standard-Pando cell types do not receive manufactured condition coefficients.

Layer 1 schema v6 is intentionally incompatible with v5 artifacts for downstream validation. Existing `step_layer1.rds` objects must be regenerated before Layer 2 so that the quantitative reaction-expression matrices are reconstructed from equal-weight mean single-cell CPM rather than the old latent-CPM estimator. Stage 1 and Stage 2 artifacts remain reusable if their existing workflow contracts are otherwise unchanged, because the new quantitative RNA estimator is reconstructed in Stage 4 from the retained Stage 1 cell-level Pando objects and the exact Stage 2 SuperCell membership.
