# RegCompass mathematical specification

This is the single canonical document for RegCompass equations and quantitative definitions. Tutorials and Rd files describe interfaces only.

## 1. Condition-comparable Pando model

For one broad cell type and target gene \(g\), Pando first applies its structural candidate rules: regulatory-domain proximity, TF motif support, and the supplied regulatory-region prior. For the human RegCompass route the default region prior is

\[
R_{prior}=R_{phastCons}\cup R_{SCREEN\ cCRE}.
\]

For exact edge \(e=(g,f,r)\), condition \(c\), define the Pando marginal correlations

\[
\rho^{peak}_{e,c}=cor_c(ATAC_r,RNA_g),\qquad
\rho^{TF}_{e,c}=cor_c(RNA_f,RNA_g).
\]

The same correlations are computed once more using all eligible-condition cells pooled within the broad cell type. The implementation uses strict comparisons. With configured thresholds \(\tau_{peak}\) and \(\tau_{TF}\),

\[
L_{e,c}=\mathbf 1\{|\rho^{peak}_{e,c}|>\tau_{peak}\}
\mathbf 1\{|\rho^{TF}_{e,c}|>\tau_{TF}\},
\]

\[
G_e=\mathbf 1\{|\rho^{peak}_{e,global}|>\tau_{peak}\}
\mathbf 1\{|\rho^{TF}_{e,global}|>\tau_{TF}\}.
\]

Both `tf_cor` and `peak_cor` default to 0.05. They are adjustable parameters, not fixed constants; the same user-supplied values are used in pooled/global and condition-specific candidate discovery and are stored in the fit contract.

Let \(E_g^{global}\) be the pooled/global candidate triples and \(E_{g,c}\) the triples supported in condition \(c\). The frozen common dictionary is

\[
E_g^{\cup}=\operatorname{unique}\left(
E_g^{global}\cup\bigcup_cE_{g,c}
\right).
\]

Deduplication is on the complete `(target, TF, region)` triple. TFs, peaks, and targets are never unioned independently and recombined.

For edge \(e\), paired cell \(i\), and condition \(c\), the Pando interaction predictor remains

\[
x_{e,i,c}=T_{e,i}A_{e,i},
\]

where \(T\) is TF expression and \(A\) is peak accessibility. Every condition uses the same ordered columns \(e\in E_g^{\cup}\). Predictors share the equal-condition RMS within-condition scaling convention

\[
s_e^2=\frac{1}{K}\sum_{c=1}^{K}
\frac{1}{n_c}\sum_{i\in c}(x_{e,i,c}-\bar x_{e,c})^2,
\]

and coefficients are back-transformed to original TF×ATAC units.

For one target and \(K\) conditions the no-fusion ridge objective is

\[
\frac{1}{K}\sum_{c=1}^{K}
\frac{\|y_c-X_c\beta_c\|_2^2}{n_c}
+\lambda\sum_{c=1}^{K}\|\beta_c\|_2^2.
\]

One target-specific \(\lambda\) is selected by condition-stratified cross-validation and used for all condition blocks. There is **no cross-condition fusion term** and no penalty shrinking \(\beta_c\) toward \(\beta_{c'}\). Thus conditions share the dictionary coordinate system, scaling convention, and tuning parameter while retaining independently estimated coefficient blocks. Ridge makes the penalized system identifiable under severe collinearity or raw rank deficiency, although individual coefficients can remain biologically non-identifiable when predictors encode nearly the same signal.

Let \(\theta_{adj}\) denote the configurable BH threshold, default 0.05. Pando condition-specific statistical support is

\[
S_{e,c}=\mathbf 1\{estimable_{e,c}\land
P^{adj}_{e,c}<\theta_{adj}\}.
\]

Because every fitted coefficient row already corresponds to an exact edge in the frozen common dictionary, the Pando active condition edge is

\[
\boxed{A_{e,c}=S_{e,c}}.
\]

The marginal-correlation quantities \(G_e\) and \(L_{e,c}\) have one role: common-dictionary admission and provenance. They are **not** applied again as a post-fit activity gate. Consequently, if an edge enters the dictionary through pooled/global support or through condition \(A\), condition \(B\) can still call that edge active when \(\widehat\beta_{e,B}\) is estimable and BH-supported even if \(G_e=0\) and \(L_{e,B}=0\). This is intentional because the marginal correlations and the coefficient of the joint TF×ATAC ridge model are different statistics; reapplying the marginal-correlation gate would reintroduce the small-condition false-negative problem that the common dictionary is designed to mitigate.

Pando exports

\[
penalty\_effect_{e,c}=A_{e,c}\widehat\beta_{e,c}.
\]

RegCompass does not reselect Pando edges. It adds target-model validity and a hard selected-\(\lambda\) full-data \(R^2\) requirement before downstream use. Let \(\tau_{R^2}\) denote `RegCompassR.target_rsq_threshold` (default 0.05). RegCompass penalty eligibility is

\[
H_{e,c}=A_{e,c}
\mathbf 1\{fit\_status_{g,c}=\text{"ok"}\}
\mathbf 1\{R^2_{g,c}\ge\tau_{R^2}\}.
\]

The complete common-dictionary coefficient table remains stored for diagnostics and direct condition contrasts even when \(H_{e,c}=0\). Direct differential inference uses \(\Delta\beta_e=\beta_{e,A}-\beta_{e,B}\); significance in one condition and nonsignificance in another is not used as a substitute for this contrast.

The target weight used by the downstream penalty is binary. A target-condition pair with at least one penalty-eligible active edge has

\[
q_{g,c}=1.
\]

An evaluated target with no penalty-eligible edge has \(q=0\); a target without a valid finite target fit is unavailable for regulatory projection. The selected-\(\lambda\) full-data \(R^2\) therefore acts **only as the target hard gate** above: it does not continuously rescale \(q\), the regulatory projection, or the COMPASS-like penalty. Out-of-fold \(R^2\) is not part of the current canonical Pando→RegCompass handoff. The same binary target rule is used by the standard-Pando route after its full-data \(R^2\) gate.

Because correlation screening, ridge estimation, and CV tuning use the observed data, ridge-Wald and contrast P values are approximate and conditional on the selected dictionary and tuning procedure; they are not exact selective-inference P values.

## 2. Metacell regulatory projection

For metacell \(u\), RegCompass averages TF expression and peak accessibility separately over the exact member cells before multiplying them. For edge \(e\),

\[
\overline T_{e,u}=\frac{1}{|M_u|}\sum_{i\in M_u}T_{e,i},\qquad
\overline A_{e,u}=\frac{1}{|M_u|}\sum_{i\in M_u}A_{e,i}.
\]

Using only penalty-eligible active edges, the target regulatory score is

\[
G_{g,u}=\sum_{e\in E_g^{\cup}}
H_{e,c(u)}\widehat\beta_{e,g,c(u)}\,\overline T_{e,u}\,\overline A_{e,u}.
\]

This is **\(\beta\times mean(TF)\times mean(ATAC)\)** for each active edge, followed by target-level summation; it is not the mean of cell-wise TF×ATAC products.

For target \(g\) within cell type \(t\), define a robust calibration scale

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},
10^{-6}
\right).
\]

The bounded regulatory modifier is

\[
R_{g,u}=q_{g,c(u)}\tanh\left(\frac{G_{g,u}}{\sigma_{g,t(u)}}\right),
\qquad -1\le R_{g,u}\le1,
\]

with \(q=1\) for a valid active target. Thus model-fit \(R^2\) values never attenuate the modifier after passing the hard target gate. When regulatory evidence is unavailable or rejected, \(R_{g,u}=0\).

## 3. Quantitative RNA input to the COMPASS-like penalty

Let \(Y_{g,i}\) be the raw RNA count of gene \(g\) in original cell \(i\), and

\[
L_i=\sum_hY_{h,i}
\]

its complete RNA library size. RegCompass first computes linear per-cell CPM,

\[
x_{g,i}=10^6\frac{Y_{g,i}}{L_i},
\]

then takes an equal-weight mean across the final SuperCell membership,

\[
X^{RNA}_{g,u}=\frac{1}{|M_u|}\sum_{i\in M_u}x_{g,i}.
\]

Thus the quantitative LP path uses `mean(single-cell CPM)` rather than `CPM(sum counts)`, and it does not use the empirical-Bayes latent CPM estimator.

The regulatory modifier acts multiplicatively,

\[
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}}.
\]

Because \(-1\le R\le1\), the regulatory multiplier is bounded between \(1/2\) and \(2\).

## 4. GPR aggregation and quantitative reaction cost

For reaction \(r\), let \(A_{r,j}\) denote one AND branch of its Boolean GPR. With the default `gpr_and_method = "min"`,

\[
Q^{quant}_{r,j,u}=\min_{g\in A_{r,j}}X^{MO}_{g,u}.
\]

`"mean"` or `"median"` replace the AND operator when explicitly selected. Isozyme/OR branches are additive,

\[
E^{quant}_{r,u}=\sum_jQ^{quant}_{r,j,u}.
\]

The reaction coefficient supplied to the COMPASS-like LP objective is

\[
p_{r,u}=\frac{1}{1+\log_2(1+\max(E^{quant}_{r,u},0))}.
\]

Missing quantitative expression and structural-only reaction roles use cost \(p=1\).

## 5. Separate bounded structural-support path

CORDA2 structural evidence uses a bounded representation distinct from the quantitative LP input. Let \(\widehat{CPM}^{struct}_{g,u}\) denote the structural-only latent metacell CPM estimate and

\[
L^{struct}_{g,u}=\ln(1+\widehat{CPM}^{struct}_{g,u}).
\]

With half-saturation parameter \(h\),

\[
C^{RNA}_{g,u}=\frac{L^{struct}_{g,u}}{L^{struct}_{g,u}+h}.
\]

Regulatory evidence modifies structural-support odds,

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{R_{g,u}}}.
\]

The same Boolean GPR topology is then applied to obtain bounded reaction structural support. This bounded quantity is used for structural-confidence classification and is **not** substituted for \(E^{quant}\) in the LP penalty.

## 6. Medium constraints

A medium scenario changes exchange bounds by intersection with the parent GEM. For reaction \(r\),

\[
l'_r\ge l_r,\qquad u'_r\le u_r.
\]

Therefore medium application cannot create a reaction direction that was blocked in the original GEM.

## 7. CORDA2 structural reconstruction

For cell type \(t\), condition-specific reaction catalogues are unioned within that cell type before structural reconstruction. Bounded evidence is mapped to CORDA2 confidence classes \(HC,MC,NC,OT\), and one medium-constrained parent GEM is reconstructed for each cell-type × medium combination.

The current canonical implementation has two additional structural rules:

1. every input core reaction is an immutable structural requirement and cannot be deleted by CORDA2 finalization;
2. after CORDA2 finishes, no second parent/final closure LP pass is run.

Internally reversible reactions are evaluated directionally by the CORDA2 state machine. During final merge, a reaction is retained when either allowed split direction is selected; all core reactions are retained unconditionally. Retained reactions recover the corresponding medium-constrained parent reaction bounds. Consequently the final CORDA2 GEM is the model passed directly to downstream COMPASS-like scoring.

FASTCORE is an explicit supplementary completion route. `model_mode = "full_gem"` bypasses context-specific reconstruction and scores the complete medium-constrained GEM.

## 8. Directional feasibility and COMPASS-like LP scoring

For final cell-type/medium model \(\mathcal G_{t,m}\), target reaction \(r\), and direction \(d\in\{+1,-1\}\), the scoring stage first computes

\[
v^{max}_{r,d,t,m}=\max_v d\,v_r
\]

subject to

\[
S_{t,m}v=0,\qquad l_{t,m}\le v\le u_{t,m}.
\]

This is the single post-reconstruction directional feasibility calculation used for scoring; it is not a second CORDA2 closure stage.

For metacell \(u\), RegCompass then solves

\[
P^*_{r,d,u,m}=\min_{v,z}\sum_jp_{j,u}z_j
\]

subject to

\[
S_{t,m}v=0,
\qquad -z_j\le v_j\le z_j,
\qquad d\,v_r\ge\omega v^{max}_{r,d,t,m}.
\]

The normalized penalty is

\[
\widetilde P_{r,d,u,m}=\frac{P^*_{r,d,u,m}}
{\omega v^{max}_{r,d,t,m}}.
\]

Lower normalized penalty indicates stronger model-constrained support for the target direction. It is not measured metabolic flux.

The comparison/ranking score is

\[
S_{r,d,u,m}=-\log(\widetilde P_{r,d,u,m}+\epsilon).
\]

Larger score indicates stronger model support.

## 9. RNA-only interpretation control

The RNA-only control reuses the exact same structural GEM, medium, bounds, target directions, and \(v^{max}\). Only the quantitative objective coefficient differs:

\[
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}}
\quad\text{versus}\quad
X^{RNA}_{g,u}.
\]

Therefore the multiome-versus-RNA-only penalty difference isolates the effect of the regulatory modifier conditional on the already constructed structural model; it is not a separately reconstructed RNA-only GEM.

## 10. Condition statistics

Condition comparisons are performed within a fixed cell type, reaction, direction, and medium. Pairwise comparisons use Wilcoxon rank-sum tests; analyses with at least three conditions can additionally use a Kruskal–Wallis omnibus test. Multiple-testing correction is applied according to the requested adjustment scope.

Metacells are within-dataset statistical units. These tests are not donor/sample-level biological-replicate inference unless the study design supplies an appropriate independent replicate level outside this metacell test.

## 11. Artifact compatibility

The quantitative LP path and bounded structural path are deliberately separate. Layer 1 schema-v6 artifacts expose both quantitative reaction expression and bounded structural support. Older Layer 1 artifacts that lack the quantitative matrices must be regenerated before canonical Layer 2 scoring.