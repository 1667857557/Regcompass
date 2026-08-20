# RegCompass mathematical specification

This document is the canonical quantitative specification. Tutorials and Rd files describe interfaces; the equations below define the production calculations.

## 1. Condition-comparable Pando E★-JSE model

For one broad cell type and target gene \(g\), Pando first applies its structural candidate rules: regulatory-domain proximity, TF motif support, and the configured peak-target and TF-target correlation screens. For the human RegCompass route the default region prior is

\[
R_{prior}=R_{phastCons}\cup R_{SCREEN\ cCRE}.
\]

For exact edge \(e=(g,f,r)\), condition \(c\), define

\[
\rho^{peak}_{e,c}=cor_c(ATAC_r,RNA_g),\qquad
\rho^{TF}_{e,c}=cor_c(RNA_f,RNA_g).
\]

The same correlations are also computed on all eligible-condition cells pooled within the broad cell type. With strict configured thresholds \(\tau_{peak}\) and \(\tau_{TF}\), pooled/global and condition-local candidate sets are unioned on the complete exact coordinate:

\[
D_T=D_{global}\cup D_1\cup\cdots\cup D_K.
\]

Deduplication is on `(target, TF, region)`. TFs, peaks and targets are never unioned independently and recombined. Every condition fits the same ordered \(D_T\).

For paired cell \(i\),

\[
x_{cie}=RNA_{ci,f(e)}\,ATAC_{ci,r(e)}.
\]

Each edge uses the equal-condition within-condition RMS scale

\[
s_e=\sqrt{\frac1K\sum_{c=1}^K
\frac1{n_c}\sum_i(x_{cie}-\bar x_{c,e})^2}.
\]

This prevents the predictor scale itself from being determined by the largest condition, while the likelihood still retains the true information difference through \(X_c^TX_c\). Each condition has its own intercept:

\[
y_c=\alpha_c\mathbf 1+X_c\beta_c+\epsilon_c,
\qquad \epsilon_c\sim N(0,\sigma_g^2I).
\]

After condition-wise centering, one common target residual scale defines

\[
Q_c=X_c^TX_c/\widehat\sigma_g^2,
\qquad h_c=X_c^Ty_c/\widehat\sigma_g^2,
\]

and

\[
Q=blockdiag(Q_1,\ldots,Q_K),\qquad
h=(h_1^T,\ldots,h_K^T)^T.
\]

There is no \(n_{total}/(Kn_c)\) condition-size equalization weight.

### 1.1 Q-orthogonal shared/deviation geometry

Let

\[
A=\mathbf1_K\otimes I_p.
\]

For each exact edge, the implementation uses a predefined reference ordering when its contrasts are identifiable. If a reference contrast is not identifiable, it constructs a maximal identifiable contrast tree rather than allowing an unidentified coefficient difference to drift. Let \(D\) be the resulting contrast operator and \(B_0\) a right inverse on the identifiable contrast subspace, \(DB_0=I\). Define

\[
R=B_0-A(A^TQA)^+A^TQB_0.
\]

Then

\[
DR=I,\qquad A^TQR=0,
\]

and the stacked coefficient has the exact decomposition

\[
\beta=A\mu+R\delta.
\]

The shared and deviation likelihood blocks are therefore Q-orthogonal. The shared component is always the joint MLE

\[
\widehat\mu=(A^TQA)^+A^Th,
\]

with no empirical-Bayes or ridge shrinkage of the common mean toward zero.

Define

\[
H=R^TQR,\qquad r=R^Th.
\]

For every identifiable contrast coordinate \(j\), its effective information \(I_j\) is obtained from the full correlated profile information, equivalently from the appropriate diagonal of \(H^+\) only after estimability is established. Production uses the fixed E★ threshold

\[
\boxed{z=0.25}
\]

and solves

\[
\widehat\delta=
\arg\min_\delta\left[
\frac12\delta^TH\delta-r^T\delta
+0.25\sum_{j\in\mathcal I}\sqrt{I_j}|\delta_j|
\right],
\]

subject to \(\delta_j=0\) for non-identifiable/zero-information coordinates. In the scalar case,

\[
\widehat\delta=sign(d)
\left(|d|-\frac{0.25}{\sqrt{I_\delta}}\right)_+.
\]

A zero-information equality is a deterministic boundary convention and is recorded as `shared_by_boundary`; it is not evidence of biological equality. An identifiable contrast shrunk exactly to zero is recorded separately as `fused_by_penalty`.

The production coefficient is

\[
\boxed{\widehat\beta^E=A\widehat\mu+R\widehat\delta}.
\]

Pando exports, in raw TF×ATAC units,

\[
estimate_{c,e}=penalty\_effect_{c,e}=\widehat\beta^E_{c,e}.
\]

The production conditional route exposes no alternative \(z\) grid, no conditional ridge-CV lambda selection and no fusion-ratio sensitivity parameter.

### 1.2 Fusion-component joint covariance and inference

Exact zero deviations define an equality graph among conditions for each edge. Connected components are the fusion components. After the E★ structure is frozen, let

\[
b=M_z\theta_z
\]

map the stacked condition-edge coefficient vector to one parameter per edge × fusion component. With condition-specific intercepts and block design \(\mathcal X=blockdiag(X_1,\ldots,X_K)\), the reduced joint design is

\[
Z_z=[I_{condition\ intercept},\;\mathcal X M_z].
\]

The inference-only selected-structure refit is

\[
\widetilde\theta_z=(Z_z^TZ_z)^+Z_z^Ty,
\]

with

\[
\widehat\sigma^2_{inf,z}
=\frac{\|y-Z_z\widetilde\theta_z\|_2^2}
{N-rank(Z_z)},
\]

and full covariance

\[
\widehat{Cov}(\widetilde\theta_z)
=\widehat\sigma^2_{inf,z}(Z_z^TZ_z)^+.
\]

For condition \(c\), edge \(e\), selection row \(L_{c,e}\),

\[
SE_{c,e}=\sqrt{L_{c,e}\widehat{Cov}(\widetilde\theta_z)L_{c,e}^T}.
\]

Thus fully or partially shared edges borrow information through the full correlated joint design rather than through diagonal-information addition. `inference_estimate` and `inference_se` are stored separately and never replace the production `estimate`/`penalty_effect`.

The Wald statistic and raw P value are

\[
Z_{c,e}=\frac{\widetilde\beta_{c,e}}{SE_{c,e}},
\qquad
p_{c,e}=2\Phi(-|Z_{c,e}|).
\]

If a reduced-model component remains non-estimable, `inference_se`, `pval`, and `padj` are `NA`; the row is not converted to `p=1` or `SE=Inf`.

### 1.3 condition × target BH and RegCompass exact-edge union

BH is performed independently for every condition-target family:

\[
q_{c,e}=BH\{p_{c,e'}:target(e')=g\}.
\]

Only finite, estimable P values enter each family, and significance is strict:

\[
S_{c,e}=1\{q_{c,e}<\theta_{adj}\},
\qquad \theta_{adj}=0.05\ \text{by default}.
\]

For exact edge \(e\), define engineering validity

\[
V_e=1\{\forall c:\ fit\_status_{c,g(e)}=ok
\land \widehat\beta^E_{c,e}\ \text{finite}\}.
\]

The RegCompass handoff is

\[
\boxed{H_e=V_e\,1\{\exists c:q_{c,e}<\theta_{adj}\}}.
\]

If \(H_e=1\), the edge is retained in **every** condition and each condition contributes its own continuous \(\widehat\beta^E_{c,e}\). If \(H_e=0\), the edge is not projected in any condition. Hence condition-specific biological effects are represented by continuous coefficient differences rather than sample-size-dependent edge presence/absence.

Pairwise contrasts are reconstructed from the same joint coefficient vector. For three conditions C/J/M,

\[
\Delta_{JM}=\Delta_{CM}-\Delta_{CJ}
\]

holds exactly. A boundary-shared contrast is not eligible for a differential-regulation claim.

The Pando-equivalent target diagnostic remains

\[
R^2_{c,g}=1-\frac{RSS_{c,g}}{TSS_{c,g}}.
\]

For conditional E★-JSE it is diagnostic only: it is not a RegCompass handoff gate and is not used to choose \(z\). Standard one-condition Pando remains a separate route and retains its own standard ridge controls and standard R² eligibility logic.

The downstream binary target availability for a condition is

\[
q_{g,c}=1
\]

when at least one admitted exact edge for target \(g\) is projected in that condition, and \(q_{g,c}=0\) when the target was validly evaluated but no exact edge was admitted. A failed/non-finite target fit remains unavailable (`NA`).

## 2. Metacell regulatory projection

For metacell \(u\), RegCompass retains the paired TF×ATAC realization from the exact same member cells and shrinks it toward the product of separate metacell means. For edge \(e\),

\[
\overline{TA}_{e,u}=\frac1{|M_u|}\sum_{i\in M_u}T_{e,i}A_{e,i},
\]

\[
\overline T_{e,u}=\frac1{|M_u|}\sum_{i\in M_u}T_{e,i},\qquad
\overline A_{e,u}=\frac1{|M_u|}\sum_{i\in M_u}A_{e,i}.
\]

The canonical exposure mixture is

\[
J_{e,u}=0.75\overline{TA}_{e,u}
+0.25\overline T_{e,u}\overline A_{e,u}.
\]

This `0.25` is the RegCompass exposure mixing weight and is mathematically independent of the E★ deviation threshold \(z=0.25\).

Using only admitted exact edges, the target regulatory score is

\[
G_{g,u}=\sum_{e\to g}H_e\widehat\beta^E_{c(u),e}J_{e,u}.
\]

The RNA and ATAC vectors used in \(\overline{TA}\) are indexed by the same paired cell IDs within the exact SuperCell membership.

For target \(g\) within cell type \(t\), define

\[
\sigma_{g,t}=\max\left(
\frac{IQR(G_{g,t})}{1.349},
MAD_{1.4826}(G_{g,t}),
\sqrt{mean(G_{g,t}^2)},10^{-6}
\right).
\]

The bounded regulatory modifier is

\[
R_{g,u}=q_{g,c(u)}\tanh\left(\frac{G_{g,u}}{\sigma_{g,t(u)}}\right),
\qquad -1\le R_{g,u}\le1.
\]

## 3. Quantitative RNA input to the COMPASS-like penalty

Let \(Y_{g,i}\) be raw RNA count and \(L_i=\sum_hY_{h,i}\) the complete cell library size. RegCompass computes single-cell linear CPM,

\[
x_{g,i}=10^6\frac{Y_{g,i}}{L_i},
\]

then the equal-weight mean within the final SuperCell membership,

\[
X^{RNA}_{g,u}=\frac1{|M_u|}\sum_{i\in M_u}x_{g,i}.
\]

The regulatory modifier acts multiplicatively,

\[
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}},
\]

so the multiplier is bounded between \(1/2\) and \(2\).

## 4. GPR aggregation and quantitative reaction cost

For reaction \(r\), let \(A_{r,j}\) denote one AND branch. With the default `gpr_and_method = "min"`,

\[
Q^{quant}_{r,j,u}=\min_{g\in A_{r,j}}X^{MO}_{g,u}.
\]

`"mean"` or `"median"` replace the AND operator when explicitly selected. Isozyme/OR branches are additive,

\[
E^{quant}_{r,u}=\sum_jQ^{quant}_{r,j,u}.
\]

The reaction coefficient supplied to the COMPASS-like LP objective is

\[
p_{r,u}=\frac1{1+\log_2(1+\max(E^{quant}_{r,u},0))}.
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
C^{MO}_{g,u}=\frac{C^{RNA}_{g,u}2^{R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{R_{g,u}}}.
\]

The same Boolean GPR topology is then applied to bounded reaction structural support. This structural quantity is not substituted for \(E^{quant}\) in the LP penalty.

## 6. Medium constraints

A medium scenario changes exchange bounds by intersection with the parent GEM. For reaction \(r\),

\[
l'_r\ge l_r,\qquad u'_r\le u_r.
\]

Therefore medium application cannot create a reaction direction that was blocked in the original GEM.

## 7. CORDA2 structural reconstruction

For cell type \(t\), condition-specific reaction catalogues are unioned within that cell type before structural reconstruction. Bounded evidence is mapped to CORDA2 confidence classes \(HC,MC,NC,OT\), and one medium-constrained parent GEM is reconstructed for each cell-type × medium combination.

Every input core reaction is an immutable structural requirement. After CORDA2 finishes, no second parent/final closure LP pass is run. Internally reversible reactions are evaluated directionally by the CORDA2 state machine; during final merge a reaction is retained when either allowed split direction is selected, while all core reactions are retained unconditionally. Retained reactions recover the corresponding medium-constrained parent bounds.

FASTCORE is an explicit supplementary completion route. `model_mode = "full_gem"` bypasses context-specific reconstruction and scores the complete medium-constrained GEM.

## 8. Directional feasibility and COMPASS-like LP scoring

For final cell-type/medium model \(\mathcal G_{t,m}\), target reaction \(r\), and direction \(d\in\{+1,-1\}\),

\[
v^{max}_{r,d,t,m}=\max_v d\,v_r
\]

subject to

\[
S_{t,m}v=0,\qquad l_{t,m}\le v\le u_{t,m}.
\]

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
{\omega v^{max}_{r,d,t,m}},
\]

and the comparison score is

\[
S_{r,d,u,m}=-\log(\widetilde P_{r,d,u,m}+\epsilon).
\]

Lower normalized penalty and higher score indicate stronger model support; neither is measured flux.

## 9. RNA-only interpretation control

The RNA-only control reuses the same structural GEM, medium, bounds, target directions, and \(v^{max}\). Only the quantitative objective coefficient differs:

\[
X^{MO}_{g,u}=X^{RNA}_{g,u}2^{R_{g,u}}
\quad\text{versus}\quad X^{RNA}_{g,u}.
\]

Thus the multiome-versus-RNA-only penalty difference isolates the regulatory modifier conditional on the same structural model.

## 10. Condition statistics

Condition comparisons are performed within a fixed cell type, reaction, direction, and medium. Pairwise comparisons use Wilcoxon rank-sum tests; analyses with at least three conditions can additionally use a Kruskal–Wallis omnibus test. Multiple-testing correction follows the requested adjustment scope.

Metacells are within-dataset statistical units; these tests are not donor/sample-level biological-replicate inference unless an appropriate independent replicate level is supplied by the study design.

## 11. Artifact compatibility

The quantitative LP path and bounded structural path are deliberately separate. Layer 1 artifacts expose both quantitative reaction expression and bounded structural support. Artifacts generated with the old conditional no-fusion ridge/local-Wald topology, an R² hard gate on conditional handoff, or the previous pure product-of-means Pando projection must be regenerated before comparison with the E★-JSE z=0.25 production workflow.
