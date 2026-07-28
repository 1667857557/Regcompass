# RegCompassR 1.8.10 workflow

## Canonical data flow

```text
all conditions within one cell type
→ one validated GREAT-domain Pando TF–peak–GEM-gene candidate universe
→ condition-aware observability filter
→ direct condition-specific theta elastic net
→ derived cross-condition backbone and zero-sum deviations
→ sample/donor cluster bootstrap within condition when sample_col is valid
→ explicit warning and cell-bootstrap fallback otherwise
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

The canonical workflow requires condition and cell-type metadata. An optional biological sample column is used only for Stage 1 bootstrap; Stage 2 remains condition-only.

## 1. Shared Pando candidate background

For cell type \(m\), Pando constructs one condition-agnostic structural set

\[
\mathcal U_m^{struct}=\{e=(t,p,g)\},
\]

where \(t\) is a measured TF, \(p\) is an exact measured ATAC feature supported by a regulatory region and motif, and \(g\) is a GEM GPR target gene. The canonical `peak_to_gene_method = "GREAT"` creates broad basal-plus-extension domains without using target-expression correlation to admit candidates. Every condition uses the same candidate dictionary and design fingerprint.

RegCompass retains candidates that have an observable `TF RNA × peak ATAC` predictor and target RNA in at least one condition. This shared model dictionary receives `model_edge_universe_id`.

## 2. Direct condition-specific sparse coefficients

For edge \(e=(t,p,g)\) and cell \(u\),

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target and predictors are centred within condition and each edge receives one scale shared across conditions:

\[
y^\circ_u=\sum_e \widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

The elastic-net objective is

\[
\min_\Theta
\frac12\sum_u w_u
\left(y^\circ_u-\sum_e\widetilde x_{e,u}\theta_{e,c(u)}\right)^2
+\lambda\alpha p_\theta\sum_{e,c}|\theta_{e,c}|
+\frac{\lambda(1-\alpha)p_\theta}{2}\sum_{e,c}\theta_{e,c}^2,
\]

with \(w_u\propto1/n_{c(u)}\). The reported backbone and deviations are derived:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e,
\qquad
\sum_c\delta_{e,c}=0.
\]

## 3. Leakage-resistant CV and sample-aware bootstrap

Five folds are assigned separately inside every condition. Training-fold condition means and edge scales are applied to validation cells. Cross-validation remains cell-level; active targets require positive out-of-fold \(R^2\).

For condition \(c\), let \(D_c\) denote the observed sample IDs. With a valid `sample_col`, each bootstrap draws \(|D_c|\) sample IDs with replacement and includes all cells from each selected sample. This preserves donor clusters, so replicate cell counts can vary.

When `sample_col` is omitted or the named column does not exist, RegCompass prints the exact reason and samples the original number of cells with replacement inside each condition. Existing incomplete sample IDs are rejected.

For successful replicates,

\[
\widehat\Pi_{e,c}=\frac1{B_s}
\sum_b I(|\widehat\theta^{(b)}_{e,c}|>\varepsilon)
\]

and

\[
\rho_{e,c}=\left|
\frac{\sum_b I^{(b)}_{e,c}\operatorname{sign}(\widehat\theta^{(b)}_{e,c})}
{\sum_b I^{(b)}_{e,c}}
\right|.
\]

Bootstrap provenance is recorded before Stage 3 reads active edges. Sample-cluster bootstrap measures reproducibility across observed samples; it is not a treatment-effect test.

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

All required AND subunits must be present; alternative isoenzymes remain OR branches. Positive and negative active edges both establish target membership, so `core` means regulatory anchoring rather than proven reaction activation.

## 5. ATAC projection and reaction evidence

Stable condition coefficients are projected onto metacell ATAC using a shared TF reference and edge scale. The bounded accessibility modifier updates RNA support on the log-odds scale. GPR AND branches use `min` by default and isozyme OR branches are additive. Reaction support becomes a monotonically decreasing LP penalty.

## 6. Shared final metabolic model

Condition-specific biological reaction catalogues are merged before model construction. For each medium, one union GEM is built and one global add-only FASTCORE completion is performed. Every condition and metacell uses identical \(S\), lower bounds, and upper bounds within that medium. Only evidence-derived penalties vary.

## 7. Execution timing

Each stage prints elapsed time and final status in R after its final artifact is committed. Timing is not persisted in stage objects, the compact result, `step_timing.tsv`, or `00_execution_timing.tsv`.
