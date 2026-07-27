# RegCompassR 1.8.9 workflow

## Canonical data flow

```text
all conditions within one cell type
→ one validated Pando TF–peak–GEM-gene candidate universe
→ condition-aware observability filter
→ direct condition-specific theta elastic net
→ derived cross-condition backbone and zero-sum deviations
→ full-size condition-stratified bootstrap stability
→ condition-specific metabolic target genes
→ complete-GPR condition core reactions
→ one ordered subsystem/cross-reference expansion pass
→ merged biological reaction catalogue
→ RNA support modified by ATAC regulatory state
→ COMPASS-compatible GPR aggregation
→ one medium-specific union GEM
→ one global FASTCORE completion
→ directional two-step LP scoring
```

The canonical workflow uses condition and cell-type metadata only. It does not
accept or interpret a biological-sample column.

## 1. Shared Pando candidate background

For cell type \(m\), Pando constructs one condition-agnostic structural set

\[
\mathcal U_m^{struct}=\{e=(t,p,g)\},
\]

where \(t\) is a measured TF, \(p\) is an exact measured ATAC feature supported
by a regulatory region and motif, and \(g\) is a GEM GPR target gene present in
the RNA assay. Every condition uses the same candidate dictionary and design
fingerprint.

RegCompass then retains candidates that have an observable
`TF RNA × peak ATAC` predictor and observable target RNA in at least one
condition. This second shared dictionary receives `model_edge_universe_id`.
No target correlation, fitted coefficient or P value is used before fitting.

## 2. Direct condition-specific sparse coefficients

For edge \(e=(t,p,g)\) and cell \(u\),

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target and predictors are centred within condition and each edge receives one
scale shared across conditions. The fitted coefficient is directly
condition-specific:

\[
y^\circ_u=\sum_e \widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

The elastic-net objective is

\[
\min_\Theta
\frac12\sum_u w_u
\left(y^\circ_u-\sum_e\widetilde x_{e,u}\theta_{e,c(u)}\right)^2
+\lambda\alpha p_\theta\sum_{e,c}|\theta_{e,c}|
+\frac{\lambda(1-\alpha)p_\theta}{2}
\sum_{e,c}\theta_{e,c}^2,
\]

with \(w_u\propto1/n_{c(u)}\). The L1 penalty therefore acts directly on
\(	heta_{e,c}\), so an edge can be exactly zero in one condition while remaining
non-zero in another.

The reported shared backbone and deviations are derived summaries:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

`global_penalty_factor` and `deviation_penalty_factor` remain compatibility
aliases for the one common \(p_\theta\) and must be equal. The old statement that
condition deviations receive stronger default shrinkage is no longer valid.

## 3. Leakage-resistant CV and bootstrap sub-GRNs

Five folds are assigned separately inside every condition. Training-fold
condition means and edge scales are applied to validation cells. The selected
lambda uses `lambda.1se`; active targets must have strictly positive out-of-fold
\(R^2\).

Each bootstrap replicate resamples the original number of cells inside every
condition and refits at the selected lambda. For successful replicates,

\[
\widehat\Pi_{e,c}=\frac1{B_s}
\sum_b I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon)
\]

and

\[
\rho_{e,c}=
\left|
\frac{\sum_b I^{(b)}_{e,c}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_b I^{(b)}_{e,c}}
\right|.
\]

An active condition edge must pass bootstrap completion, selection-frequency,
sign-stability, effect-size and CV gates. Because sparsity is applied directly
to \(	heta\), binary condition sub-GRN differences are now part of the fitted
model rather than being produced only by post-fit cancellation and thresholds.

## 4. Condition genes and complete-GPR cores

The condition target set is

\[
G_{m,c}=\{g:\exists e=(t,p,g)\text{ active in }c\}.
\]

For reaction \(r\) with alternative GPR branches \(B_{r,k}\),

\[
Core_{r,m,c}=1
\iff
\exists k:B_{r,k}\subseteq G_{m,c}.
\]

All required AND subunits must be present; alternative isoenzymes remain OR
branches. Positive and negative active regulatory edges both establish target
membership, so `core` means regulatory anchoring rather than proven reaction
activation.

## 5. ATAC projection and reaction evidence

Stable condition coefficients are projected onto metacell ATAC using a shared TF
reference and edge scale. The bounded accessibility modifier updates RNA support
on the log-odds scale. GPR AND branches use `min` by default and isozyme OR
branches are additive. Reaction support becomes a monotonically decreasing LP
penalty.

## 6. Shared final metabolic model

Condition-specific biological reaction catalogues are merged before model
construction. For each medium, one union GEM is built and one global add-only
FASTCORE completion is performed. Every condition and metacell therefore uses
identical \(S\), lower bounds and upper bounds within that medium. Only
evidence-derived penalties vary.

This design supports same-reaction comparisons across conditions, but the
resulting scores remain evidence-compatible reaction potentials rather than
measured metabolic fluxes.
