# Shared-background direct-theta GRN mathematics and object contracts

RegCompassR 1.8.9 separates the structural regulatory hypothesis, the shared
observable model dictionary, and the condition-specific sparse network.

## 1. Shared structural background

For each cell type \(m\), all conditions are combined only to construct one
Pando structural dictionary

\[
\mathcal U_m^{struct}=\{e=(t,p,g)\}.
\]

`t` is a TF, `p` is an exact measured ATAC feature supported by a Pando region
and motif, and `g` is a GEM GPR target gene. Candidate edges are deduplicated by
`(TF, ATAC feature, target)` and receive `edge_universe_id`.

Pooled Pando detection thresholds are zero. RegCompass subsequently applies a
condition-aware detection-only filter to form
\(\mathcal U_m^{model}\), which receives `model_edge_universe_id`. Every
condition of the cell type uses exactly this model dictionary.

## 2. Predictor, centring and shared scale

For cell \(u\) and edge \(e=(t,p,g)\),

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target and predictor values are centred within condition:

\[
y^\circ_{g,u}=y_{g,u}-\bar y_{g,c(u)},
\qquad
x^\circ_{e,u}=x_{e,u}-\bar x_{e,c(u)}.
\]

The edge scale is equal-condition weighted and shared across conditions:

\[
s_{e,m}=\sqrt{\frac1C\sum_c\frac1{n_c}
\sum_{u\in c}(x^\circ_{e,u})^2},
\qquad
\widetilde x_{e,u}=x^\circ_{e,u}/s_{e,m}.
\]

## 3. Direct sparse condition coefficients

RegCompass fits one coefficient per edge and condition:

\[
y^\circ_u=\sum_e\widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

With condition-balanced weights \(w_u\propto1/n_{c(u)}\), the objective is

\[
\begin{aligned}
\min_\Theta\quad
&\frac12\sum_uw_u
\left(y^\circ_u-\sum_e\widetilde x_{e,u}\theta_{e,c(u)}\right)^2\\
&+\lambda\alpha p_\theta\sum_{e,c}|\theta_{e,c}|\\
&+\frac{\lambda(1-\alpha)p_\theta}{2}
\sum_{e,c}\theta_{e,c}^2.
\end{aligned}
\]

The L1 non-smooth points are therefore located at
\(\theta_{e,c}=0\). A regulatory edge can be selected in one condition and
exactly absent in another without requiring cancellation between a global and a
deviation coordinate.

The global backbone is reported as the cross-condition mean:

\[
\beta_e=\frac1C\sum_c\theta_{e,c}.
\]

Condition deviations are derived as

\[
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

Thus `global_estimate` remains a comparable shared summary, while sparsity is
applied where the biological topology is interpreted: the condition edge.

## 4. Parameter compatibility

`global_penalty_factor` and `deviation_penalty_factor` are retained as
compatibility aliases for the common direct-theta penalty \(p_\theta\). They
must be equal. The previous design that used a stronger separate deviation
penalty is retired because it did not directly induce condition-specific edge
absence.

Canonical values are

\[
\alpha=0.5,
\qquad
p_\theta=1.
\]

## 5. Leakage-resistant cross-validation

Five folds are stratified within condition. Condition means and edge scales are
estimated on training cells and applied to validation cells. One common lambda
path and one `lambda.1se` selection rule are used for the joint direct-theta
model.

The reported

\[
R^2_{CV}=1-
\frac{\sum_uw_u(y_u^\circ-\widehat y_u^\circ)^2}
{\sum_uw_u(y_u^\circ-\bar y_w^\circ)^2}
\]

is target-model predictive reliability. Active edges require
\(R^2_{CV}>0\). It is not a donor-level significance test.

## 6. Condition-stratified bootstrap

Every replicate samples the original number of cells with replacement inside
each condition, re-centres the resampled matrices and fits the direct-theta
model at the selected lambda.

\[
\widehat\Pi_{e,c}=\frac1{B_s}\sum_b
I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon)
\]

is the frequency with which the condition coefficient itself is non-zero.
Sign stability is

\[
\rho_{e,c}=\left|
\frac{\sum_bI^{(b)}_{e,c}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_bI^{(b)}_{e,c}}
\right|.
\]

The downstream coefficient is

\[
\widetilde\theta_{e,c}=
\widehat\theta_{e,c}\widehat\Pi_{e,c}\rho_{e,c}.
\]

## 7. Condition genes and complete-GPR cores

The condition target set is

\[
G_{m,c}=\{g:\exists e=(t,p,g)\text{ active in }c\}.
\]

For reaction \(r\) with GPR branches \(B_{r,k}\),

\[
Core_{r,m,c}=1
\iff
\exists k:B_{r,k}\subseteq G_{m,c}.
\]

Required complex subunits remain AND relationships and alternative isoenzymes
remain OR relationships.

## 8. ATAC-only regulatory projection

For stable edge coefficient \(\widetilde\theta_{tpg,c}\), the ATAC projection
weight is

\[
\omega_{tpg,c}=\widetilde\theta_{tpg,c}
\frac{\bar T_{t,m}}{s_{tpg,m}}.
\]

TFs sharing one measured peak and target are signed-summed. A bounded ATAC
deviation modifier updates RNA support on the log-odds scale before GPR
aggregation.

## 9. Shared final metabolic structure

Condition core and biological reaction catalogues may differ, but they are
merged before model construction. For each medium, one union GEM undergoes one
global FASTCORE completion and is reused for every condition and metacell:

\[
S_A=S_B=S^{union},
\qquad
lb_A=lb_B,
\qquad
ub_A=ub_B.
\]

Only evidence-derived reaction penalties vary. This separates condition
regulatory evidence from structural model differences.

## 10. Output contract

Stage 1 retains:

```text
tf_peak_gene_candidates
tf_peak_gene_global
tf_peak_gene_condition_all
tf_peak_gene_significant
condition_target_genes
target_model_diagnostics
stability_diagnostics
group_status
celltype_fit_status
```

Condition tables include `coefficient_parameterization =
"direct_condition_theta"`, `theta_penalty_factor`, direct
`effective_estimate`, derived `global_estimate` and `condition_deviation`, and
bootstrap/CV diagnostics.
