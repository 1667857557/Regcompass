# Shared-background direct-theta GRN mathematics and object contracts

RegCompassR 1.8.10 separates the structural regulatory hypothesis, the shared observable model dictionary, the condition-specific sparse network, and the bootstrap resampling unit.

## 1. Shared structural background

For each cell type \(m\), all conditions are combined to construct one Pando structural dictionary

\[
\mathcal U_m^{struct}=\{e=(t,p,g)\}.
\]

`t` is a TF, `p` is an exact measured ATAC feature supported by a Pando region and motif, and `g` is a GEM GPR target gene. Candidate edges are deduplicated by `(TF, ATAC feature, target)` and receive `edge_universe_id`.

Pooled Pando detection thresholds are zero. RegCompass subsequently applies a condition-aware detection-only filter to form \(\mathcal U_m^{model}\), which receives `model_edge_universe_id`. Every condition of the cell type uses exactly this model dictionary.

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

The L1 non-smooth points are located at \(\theta_{e,c}=0\). A regulatory edge can be selected in one condition and exactly absent in another without cancellation between a global and deviation coordinate.

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

## 4. Parameter compatibility

`global_penalty_factor` and `deviation_penalty_factor` are compatibility aliases for the common direct-theta penalty \(p_\theta\). They must be equal. Canonical values are

\[
\alpha=0.5,
\qquad
p_\theta=1.
\]

## 5. Leakage-resistant cross-validation

Five folds are stratified within condition. Condition means and edge scales are estimated on training cells and applied to validation cells. One common lambda path and one `lambda.1se` selection rule are used for the direct-theta model.

The reported

\[
R^2_{CV}=1-
\frac{\sum_uw_u(y_u^\circ-\widehat y_u^\circ)^2}
{\sum_uw_u(y_u^\circ-\bar y_w^\circ)^2}
\]

is target-model predictive reliability. Active edges require \(R^2_{CV}>0\). Cross-validation remains cell-level and is not a donor-level treatment-effect test.

## 6. Condition-stratified bootstrap

Let \(D_c\) be the observed sample/donor identifiers in condition \(c\). When a valid `sample_col` is supplied, bootstrap replicate \(b\) draws a multiset

\[
D_c^{(b)}=(d_{c,1}^{(b)},\ldots,d_{c,|D_c|}^{(b)}),
\qquad
d_{c,j}^{(b)}\overset{iid}{\sim}D_c.
\]

Every occurrence of selected donor \(d\) contributes all cells belonging to that donor and condition:

\[
I_c^{(b)}=
\mathop{\Vert}_{d\in D_c^{(b)}}
\{u:c(u)=c,\ d(u)=d\},
\]

where \(\Vert\) denotes concatenation with multiplicity. This is a cluster bootstrap. Donor cluster sizes are preserved, so the number of cells in one replicate need not equal the original \(n_c\).

If `sample_col` is omitted or the named column does not exist, RegCompass prints the exact fallback reason and samples \(n_c\) cells with replacement inside each condition. An existing sample column containing missing or empty identifiers is an error.

Both bootstrap modes re-centre the resampled target and predictors within condition and fit the direct-theta model at the full-data selected lambda and edge scale.

Selection frequency is

\[
\widehat\Pi_{e,c}=\frac1{B_s}\sum_b
I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon).
\]

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

Sample-cluster bootstrap estimates reproducibility with respect to the observed biological sample clusters. It does not by itself estimate a condition treatment effect, and reliability is weak when a condition contains very few samples. Fallback cell bootstrap estimates cell-resampling reproducibility only.

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

Required complex subunits remain AND relationships and alternative isoenzymes remain OR relationships.

## 8. ATAC-only regulatory projection

For stable edge coefficient \(\widetilde\theta_{tpg,c}\), the ATAC projection weight is

\[
\omega_{tpg,c}=\widetilde\theta_{tpg,c}
\frac{\bar T_{t,m}}{s_{tpg,m}}.
\]

TFs sharing one measured peak and target are signed-summed. A bounded ATAC deviation modifier updates RNA support on the log-odds scale before GPR aggregation.

## 9. Shared final metabolic structure

Condition core and biological reaction catalogues may differ, but they are merged before model construction. For each medium, one union GEM undergoes one global FASTCORE completion and is reused for every condition and metacell:

\[
S_A=S_B=S^{union},
\qquad
lb_A=lb_B,
\qquad
ub_A=ub_B.
\]

Only evidence-derived reaction penalties vary. The bootstrap resampling unit can alter selected regulatory edges and condition cores, but it does not create condition-specific GEM structures.

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
bootstrap_policy
```

Condition and stability tables include:

```text
bootstrap_method
bootstrap_resampling_unit
bootstrap_sample_col
n_bootstrap_samples_total
min_bootstrap_samples_per_condition
bootstrap_fallback_reason
selection_frequency
sign_stability
```

The compact final result records the selected `sample_col`, bootstrap resampling unit, fallback reason, and Stage 1 bootstrap policy in result parameters/provenance. Execution time is printed in R and is not retained in result objects or timing files.
