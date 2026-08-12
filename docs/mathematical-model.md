# RegCompass mathematical specification

This is the single canonical document for RegCompass equations and quantitative definitions. Tutorials and Rd files describe interfaces only.

## 1. Condition-comparable Pando model

For one broad cell type and target gene \(g\), let \(E_g^{pool}\) be the candidate TF–peak–target triples discovered in the pooled cell type and \(E_{g,c}\) the triples discovered in condition \(c\). The preliminary common candidate dictionary is

\[
D_g^{(0)}=E_g^{pool}\cup\bigcup_c E_{g,c}.
\]

For edge \(e\), paired cell \(i\), and target \(g\), the Pando interaction predictor is

\[
x_{e,i}=T_{e,i}A_{e,i},
\]

where \(T\) is TF expression and \(A\) is peak accessibility. All retained conditions use the same ordered predictor dictionary and are fitted jointly by the multi-task ridge model. The regularized objective contains condition-balanced squared error, ordinary ridge shrinkage, and the configured fusion penalty between condition coefficient vectors. The regularization parameter \(\lambda\) is chosen by cross-validation.

### 1.1 Preliminary dictionary screening

The full candidate dictionary \(D_g^{(0)}\) is first fitted with approximate ridge-Wald inference enabled. For an estimable condition-edge coefficient,

\[
z_{e,g,c}^{screen}
=
\frac{\widehat\beta_{e,g,c}^{screen}}
{SE(\widehat\beta_{e,g,c}^{screen})},
\]

and

\[
p_{e,g,c}^{screen}
=
2\Phi\left(-\left|z_{e,g,c}^{screen}\right|\right).
\]

Within each condition, finite estimable screening P values are adjusted once by Benjamini-Hochberg. Let \(\theta_{adj}\) denote `padj_threshold`; the RegCompass default is 0.05. Condition-specific screening support is

\[
S_{e,g,c}
=
\mathbf 1\left\{
q_{e,g,c}^{screen}<\theta_{adj}
\right\}.
\]

The final shared dictionary contains an edge when at least one condition supports it:

\[
D_g^{(1)}
=
\left\{
e\in D_g^{(0)}:\sum_c S_{e,g,c}>0
\right\}.
\]

Thus `adjust_method` and `padj_threshold` define the preliminary dictionary screen and condition-specific support mask. They are not a second post-selection significance test.

### 1.2 Final common-dictionary effect refit

After \(D_g^{(1)}\) is frozen, every condition is jointly refit on exactly that same dictionary with the same multi-task ridge/fusion estimator. Coefficient inference is disabled in this final refit:

\[
\widehat\beta_{e,g,c}^{final}
=
\operatorname{MultiTaskRidgeRefit}(D_g^{(1)}),
\]

but no final coefficient \(SE\), Wald statistic, P value, or BH-adjusted P value is calculated. The preliminary `screen_pval` and `screen_padj` therefore must not be interpreted as P values for \(\widehat\beta^{final}\).

The regulatory effect passed to RegCompass is

\[
\theta_{e,g,c}
=
\begin{cases}
\widehat\beta_{e,g,c}^{final}, &
S_{e,g,c}=1,\ \widehat\beta_{e,g,c}^{final}\text{ estimable},\ 
\text{and target }fit\_status=\text{"ok"},\\
0, & \text{otherwise}.
\end{cases}
\]

This separates statistical discovery from quantitative effect estimation while keeping every condition on one final predictor space. Final condition contrasts are quantitative coefficient differences,

\[
\Delta_{e,a,b}
=
\widehat\beta_{e,a}^{final}-\widehat\beta_{e,b}^{final},
\]

without a second post-selection contrast P value.

`screen_pval` and `screen_padj` are approximate ridge-Wald screening quantities conditional on the preliminary candidate dictionary, CV-selected regularization, and fusion structure. Candidate discovery is itself data dependent, so this is not claimed to be selective-inference-corrected FDR for the complete candidate-generation procedure.

The target weight used by the downstream penalty is binary. A target-condition pair with at least one final active edge has

\[
q_{g,c}=1.
\]

Targets without a valid active edge have no regulatory contribution. Target-level \(R^2\) and out-of-fold \(R^2\) remain fit diagnostics only; they do not enter \(q\), do not rescale the regulatory projection, and do not enter the COMPASS-like penalty. The standard-Pando route retains its own one-stage adjusted-P active-edge gate and the same binary target-weight rule.

## 2. Metacell regulatory projection

For metacell \(u\), RegCompass averages TF expression and peak accessibility separately over the exact member cells before multiplying them. For edge \(e\),

\[
\overline T_{e,u}=\frac{1}{|M_u|}\sum_{i\in M_u}T_{e,i},\qquad
\overline A_{e,u}=\frac{1}{|M_u|}\sum_{i\in M_u}A_{e,i}.
\]

The target regulatory score is therefore

\[
G_{g,u}=\sum_{e\in D_g^{(1)}}
\theta_{e,g,c(u)}\,\overline T_{e,u}\,\overline A_{e,u}.
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

with \(q=1\) for a valid active target. Thus model-fit \(R^2\) values never attenuate the modifier. When regulatory evidence is unavailable, \(R_{g,u}=0\).

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
