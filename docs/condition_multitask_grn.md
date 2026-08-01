# Unified cell-type condition GRN architecture

## Scope

RegCompass no longer fits a separate Pando model for every `condition × cell type` subset. The complete normalized multiome object is passed once to `Pando::initiate_grn()` and `Pando::find_motifs()`, followed by one call to `Pando::infer_condition_grn()`.

Pando then fits each cell type across all of its conditions using one candidate TF–peak–target dictionary and one sparse-group multi-task objective. RegCompass consumes the resulting universal and condition-level Pando `Network` objects without reconstructing condition GRNs independently.

The implementation is pinned to `1667857557/Pando_regcompass@6f42c8143bec6610b001e714a51627337f6d9ba9`.

## Joint condition model

For one cell type, target gene `g`, edge `e = (TF, peak, g)`, condition `c`, and cell or aggregate unit `i`, Pando uses the interaction predictor

```text
X_eic = RNA_TF,ic × ATAC_peak,ic
```

and jointly estimates the condition coefficient matrix `B = {beta_ec}`. In schematic form, the fitted objective is

```text
sum_c w_c / 2 ||y_gc - a_gc - X_gc beta_gc||²
+ lambda (1 - alpha) / 2 ||B_g||²_F
+ lambda alpha [(1 - eta) sum_e ||B_ge,.||_2
                + eta sum_e,c |beta_gec|]
```

where `eta` is Pando's `condition_mix`. The row-group term couples the same TF–peak–target edge across conditions, while the elementwise term permits condition selection and sign changes.

This is not a block of independently fitted condition models. Conditions share the edge dictionary, scaling policy, lambda path, and the joint sparse-group penalty.

## Three coefficient views

For every cell type and TF–peak–target edge:

```text
beta_universal,e = mean_c(beta_condition,e,c)
Delta_e,c = beta_condition,e,c - beta_universal,e
```

RegCompass exports:

- `tf_peak_gene_universal`: cell-type universal coefficients.
- `tf_peak_gene_condition`: active condition-level coefficients.
- `tf_peak_gene_condition_effect`: active `Delta_e,c` coefficients.

The legacy fields `tf_peak_gene_all` and `tf_peak_gene_significant` remain internal aliases of the condition-level tables so the existing Stage 3 complete-GPR machinery receives the new unified-scale condition GRN rather than separately fitted networks.

## Analysis routing

### Condition GRN → metabolic target and core-reaction path

A condition-level edge is active when its coefficient magnitude exceeds the numerical activity threshold and its target model passes `min_model_rsq`.

For each `condition × cell type`, active condition-level edges define the supported metabolic target-gene set. A reaction becomes core only when at least one complete GPR branch is contained in that set. Shared and condition-specific edges therefore contribute to biological reaction membership according to the actual condition GRN.

### Condition effect → TF×ATAC penalty path

Condition effects do not redefine the structural reaction catalogue. They modify the evidence supplied to the reaction penalty calculation.

For metacell `u`, edge activity is

```text
A_eu = RNA_TF,u × ATAC_peak,u
```

Within each cell type, edge activity is robustly centered and scaled across all conditions and bounded by `tanh`:

```text
Z_eu = tanh[(A_eu - median_u(A_eu)) / robust_scale_e]
```

For target gene `g` in condition `c`, the signed modifier is

```text
R_gu = reliability_gc ×
       sum_e q_ec sign(Delta_e,c) Z_eu

q_ec = |Delta_e,c| / sum_e |Delta_e,c|
```

`reliability_gc` is derived from the finite condition-network target `R²`. The modifier is bounded to `[-1, 1]`.

The target metabolic-gene RNA support remains the baseline:

```text
C_multiome = C_RNA × 2^(alpha R) /
             [1 - C_RNA + C_RNA × 2^(alpha R)]
```

This update is zero preserving and bounded. It changes the support log-odds but cannot create support when target-gene RNA support is zero.

## Why TF RNA is admitted

Pando's predictive feature is the TF-expression by peak-accessibility interaction, not peak accessibility alone. The previous RegCompass projection retained the fitted Pando coefficient but discarded the TF-expression component at the metacell projection stage.

The new projection restores the fitted predictor. RegCompass also records `tf_metabolic_target_overlap` and the modifier-level `tf_target_overlap` diagnostic. In the intended use case, regulatory TF genes and metabolic target genes do not overlap, so TF RNA supplies regulatory-state information distinct from the target metabolic-gene RNA baseline.

## Statistical semantics

The condition-aware Pando solver is regularized and does not provide per-edge adjusted P values. Edge admission therefore uses coefficient activity plus target-model `R²`; `padj_threshold` and `require_padj` are retained only as recorded legacy inputs and are not used for the multi-task edge decision.

The resulting modifier is predictive regulatory evidence. It is not an independent validation measurement and is not a flux estimate. Downstream forward and reverse LP targets remain separate counterfactual penalty calculations; their combination must not be interpreted as net flux.

## Output files

Stage 1 writes:

- `pando_tf_peak_gene_universal.tsv.gz`
- `pando_tf_peak_gene_condition_all.tsv.gz`
- `pando_tf_peak_gene_condition_active.tsv.gz`
- `pando_tf_peak_gene_condition_effect_all.tsv.gz`
- `pando_tf_peak_gene_condition_effect_active.tsv.gz`
- `pando_condition_network_index.tsv.gz`
- `pando_condition_fit_diagnostics.tsv.gz`, when available
- `pando_objects/condition_multitask_grn.rds`, when object saving is enabled

Layer 1 records `regulatory_mode = condition_effect_coefficient_x_TF_RNA_x_peak_ATAC` and exports the bounded gene-level regulatory modifier used for the reaction-capacity calculation.
