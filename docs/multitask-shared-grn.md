# Shared-backbone multitask GRN model

RegCompassR 1.8.8 uses one condition-comparable regulatory design per cell type and separates structural model comparability from condition-specific regulatory and metabolic evidence.

## 1. Shared Pando candidate universe

For cell type `m`, all conditions are pooled only for structural candidate construction:

```text
D_m = union over conditions D_m,c
```

Pando returns:

```text
U_m = {e = (TF_t, peak_p, target_g)}
```

An edge enters the default structural universe when:

1. `peak_p` lies in the candidate regulatory domain of `target_g`;
2. the peak contains a motif assigned to `TF_t`;
3. TF RNA, peak ATAC and target RNA satisfy minimum detection filters;
4. the source RNA and ATAC feature IDs satisfy the exported matrix contract.

No pooled target-RNA p-value is required in structural mode. This is necessary because

```text
theta[e,A] > 0 and theta[e,B] < 0
```

can produce a pooled correlation near zero even when both condition effects are strong.

The optional `union_within_group_correlation` screen computes maximum absolute within-condition peak–target and TF–target correlations. The resulting candidate universe remains common to every condition.

## 2. Edge feature

For cell or metacell `u`:

```text
x[e,u] = T[t,u] * A[p,u]
```

where `T` is normalized TF RNA and `A` is normalized peak accessibility. The product is formed before scaling so the feature retains the Pando interpretation that TF availability and binding-site accessibility jointly gate the edge.

For every edge, the common cell-type scale is

```text
s[e,m] = sqrt(mean over conditions of Var_c(x[e,.])))
```

and the fitted feature is

```text
x_tilde[e,u] = x[e,u] / max(s[e,m], epsilon)
```

A separate scale is not estimated for each condition.

## 3. Global backbone and condition deviations

For target gene `g` and condition `c`:

```text
theta[e,m,c] = beta[e,m] + delta[e,m,c]
```

The stacked design includes:

- unpenalized condition intercepts, or condition-by-sample intercepts when `sample_col` is supplied;
- one global interaction block `X`;
- one condition interaction block `I(condition=c) X` per condition.

The elastic-net objective is equivalent to

```text
minimize
  sum_c 1/(2*C*n_c) * ||y_g,c - intercept_g,c - X_g,c(beta_g + eta_g,c)||^2
  + lambda * alpha * [pG ||beta_g||_1 + pD sum_c ||eta_g,c||_1]
  + lambda * (1-alpha)/2 *
      [pG ||beta_g||_2^2 + pD sum_c ||eta_g,c||_2^2]
```

with

```text
pD > pG
0 <= alpha < 1
```

The observation weights make every condition contribute the same total loss. The stronger deviation penalty shrinks weak condition effects toward the shared backbone. A non-zero ridge component gives a unique solution despite the redundant global and condition-block parameterization.

The internal fitted coefficients are converted to the reference-free reported decomposition:

```text
theta[e,c] = internal_global[e] + internal_condition[e,c]

beta[e] = mean_c(theta[e,c])

delta[e,c] = theta[e,c] - beta[e]
```

Therefore:

```text
sum_c delta[e,c] = 0
```

and changing the condition reference level does not change `theta`, `beta` or `delta`.

## 4. Cross-validation and stability selection

Cross-validation folds are stratified by condition. At the selected `lambda`, stratified subsampling estimates:

```text
Pi[e,c]  = fraction of successful resamples with theta[e,c] != 0
rho[e,c] = absolute mean coefficient sign among selected resamples
```

The continuous stability weight is

```text
r[e,c] = Pi[e,c] * rho[e,c]
```

and the exported stability-weighted edge estimate is

```text
estimate[e,c] = theta[e,c] * r[e,c]
```

An active edge must pass:

- target cross-validated R-squared;
- selection-frequency threshold;
- sign-stability threshold;
- absolute effective-effect threshold.

These are regularized-network stability diagnostics, not classical p-values. The multitask output therefore leaves `padj` as missing and records `evidence_type = "multitask_stability_selected"`.

A direction reversal is reported only when the same edge is active in at least two conditions and the active effective coefficients include both positive and negative signs.

## 5. Condition-specific regulated genes

For condition `c` and cell type `m`:

```text
G[m,c] = {g : at least one active edge (TF, peak, g) exists}
```

Both positive and negative stable edges identify a regulated target. This definition describes regulatory membership, not expression direction.

## 6. Complete-GPR condition cores

Suppose reaction `r` has GPR branches

```text
B[r,1] OR B[r,2] OR ... OR B[r,K]
```

where each branch is an AND set of required genes. The condition core indicator is

```text
Core[r,m,c] = 1
```

only when

```text
there exists k such that B[r,k] is a subset of G[m,c]
```

A partial enzyme complex is retained only as a diagnostic candidate and does not become a core reaction.

Stage 3 then applies one non-recursive annotation expansion:

```text
core subsystem
→ direct KEGG/Reactome reaction equivalence
→ direct master-Rhea reaction equivalence
```

## 7. Shared Layer 1 ATAC projection

For an active TF–peak–target edge, the ATAC projection weight is

```text
omega[t,p,g,c] = theta[t,p,g,c]
                   * stability_weight[t,p,g,c]
                   * shared_TF_reference[t,m]
                   / interaction_scale[t,p,g,m]
```

The TF reference is the mean of condition-specific mean TF expression values, so conditions contribute equally. Metacell-specific TF RNA is not reintroduced during the downstream ATAC projection.

TF edges sharing the same peak and target are signed-summed:

```text
psi[p,g,c] = sum_t omega[t,p,g,c]
```

This prevents one ATAC peak from being counted repeatedly merely because it contains motifs for several TFs.

For normalized peak deviation `D[p,u]`, define one denominator shared by all conditions:

```text
L[g,m] = max_c sum_p |psi[p,g,c]|
```

The target reliability is

```text
q[g,m] = sqrt(clip(CV_R2[g,m], 0, 1))
```

and the regulatory modifier is

```text
R[g,u,c] = q[g,m] * clip(
             sum_p psi[p,g,c] D[p,u] / (L[g,m] + epsilon),
             -1,
             1
           )
```

This shared denominator preserves relative regulatory strength between conditions. Legacy condition-Pando results retain their original per-condition normalization through a separate compatibility dispatcher.

## 8. RNA integration and GPR capacity

RNA support is bounded:

```text
C_RNA[g,u] = L[g,u] / (L[g,u] + h)
```

ATAC changes the support log-odds:

```text
C_MO[g,u,c] = C_RNA[g,u,c] * 2^(alpha * R[g,u,c]) /
              (1 - C_RNA[g,u,c] +
               C_RNA[g,u,c] * 2^(alpha * R[g,u,c]))
```

The zero-evidence identity is exact:

```text
R[g,u,c] = 0  =>  C_MO[g,u,c] = C_RNA[g,u,c]
```

For one GPR AND branch:

```text
Q[r,k,u] = min over g in B[r,k] of C_MO[g,u]
```

by default. Isozyme OR branches are additive:

```text
E[r,u] = sum_k Q[r,k,u]
```

The reaction penalty is

```text
p[r,u] = 1 / (1 + log2(1 + max(E[r,u], 0)))
```

## 9. Identical final metabolic model structure

Condition core and expanded reaction sets may differ:

```text
Core[A] != Core[B]
```

but the final catalogue is the union:

```text
Core_union = union over m,c Core[m,c]
Module_union = union over m,c MetaModule[m,c]
```

For every medium scenario, Stage 5 constructs one model:

```text
S_union, lb_union, ub_union
```

and reuses the exact matrices and bounds for every condition and metacell. Conditions differ only through their evidence-derived reaction penalties.

The directional two-step LP first computes structural maximum flux:

```text
Vmax[r,d] = maximize d * v[r]
subject to S_union v = 0 and lb_union <= v <= ub_union
```

and then computes the minimum evidence-weighted absolute flux required to maintain a fraction `omega` of `Vmax`:

```text
minimize sum_j p[j,u] z[j]
subject to
  S_union v = 0
  lb_union <= v <= ub_union
  z[j] >= v[j]
  z[j] >= -v[j]
  d * v[r] >= omega * Vmax[r,d]
```

Because the structural model is shared, reaction scores remain comparable across conditions.

## 10. Stage object contract

Stage 1 exports:

```text
tf_peak_gene_candidates
tf_peak_gene_global
tf_peak_gene_condition_all
tf_peak_gene_significant
condition_target_genes
celltype_fit_status
sample_status
design_ids
```

`tf_peak_gene_condition_all` contains:

```text
edge_id
design_id
tf
region
atac_feature_id
target
global_estimate
condition_deviation
effective_estimate
selection_frequency
sign_stability
stability_weight
estimate
active_edge
effective_direction
sign_flip_flag
cv_rsq
```

Stage 3 consumes only active condition edges to define regulated genes and complete-GPR cores. Stage 4 consumes the same active edge table for shared-scale ATAC projection. Stage 5 consumes the merged reaction catalogue and the complete Layer 1 reaction-expression matrix.

## 11. Interpretation boundary

The GRN coefficients are regularized conditional associations learned from the same multiome experiment. Condition-centred intercepts, ATAC-only downstream projection and exact RNA-only fallback reduce direct duplicate weighting, but they do not make the learned network independent causal evidence. Replicate-aware external validation remains necessary for causal claims.
