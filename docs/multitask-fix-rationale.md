# Multitask GRN corrective design note

This change resolves three correctness and maintenance problems in the canonical
multitask GRN implementation.

## One authoritative load contract

The previous package load path rebound the same internal defaults and fitter
names across several Collate-ordered files. The corrected layout is:

```text
multitask_grn.R                  numerical utilities only
multitask_grn_cv_contract.R      direct condition-theta implementation
multitask_grn_inference.R        cell-type orchestration
multitask_grn_parameter_policy.R sole defaults, validator, and final fitter
multitask_grn_output_contract.R  exported diagnostics/provenance
```

The bootstrap, ordering, and active-edge shadow-wrapper files were removed.

## Direct condition-specific sparsity

The fitted coefficients are now \(\theta_{e,c}\) themselves:

\[
y_u^\circ=\sum_e\widetilde x_{e,u}\theta_{e,c(u)}+\varepsilon_u.
\]

Elastic-net sparsity acts directly on every edge-condition coefficient. The
reported backbone and deviations are derived after fitting:

\[
\beta_e=\frac1C\sum_c\theta_{e,c},
\qquad
\delta_{e,c}=\theta_{e,c}-\beta_e.
\]

This permits an exact zero in one condition without requiring cancellation
between separately penalised latent coordinates.

## Predictor-label ordering

Candidate edges are sorted exactly once before TF RNA and peak ATAC matrices are
extracted. Matrix columns are named by `edge_id`, and shuffled-candidate
invariance is tested. This removes the previous risk that edge rows and
TF-by-peak predictor columns could be reordered independently.

## Canonical GREAT structural domains

The canonical Pando structural rule is now `peak_to_gene_method = "GREAT"` with
`extend = 1000000`. GREAT defines a broad distal structural hypothesis without
using fitted target-expression correlation. Condition-aware observability,
direct-theta elastic net, leakage-resistant CV, and bootstrap stability provide
the downstream evidence filters. `Signac` remains an explicit structural
sensitivity option.
