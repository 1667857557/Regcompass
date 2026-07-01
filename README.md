# RegCompassR

RegCompassR is a simplified Layer 1 tool for **multiome-supported, GPR-aware, sample-aware reaction capacity potential** analysis from annotated Seurat v4 RNA+ATAC objects.

The package intentionally keeps the main analysis narrow:

```text
Seurat v4 RNA+ATAC counts
→ sample × cell_type × optional condition micropools
→ RNA raw-count pool sum
→ log2(CPM + 1)
→ robust z-score by gene across pools
→ sigmoid gene score
→ sqrt promiscuity correction
→ GPR capacity:
   AND = Boltzmann minimum-biased average, tau = 0.20
   OR  = sum across isoenzyme groups
→ cell-type Q95 continuous shrinkage
→ reaction confidence = median gene confidence × observed GPR gene fraction
→ C_raw, C_rel, reaction_confidence, minimal diagnostics
```

It does **not** perform clustering, WNN construction, full Human-GEM QP, FVA, thermodynamic modeling, causal regulator discovery, or true flux inference.

## Expected input

An already annotated Seurat v4/Signac multiome object with:

- RNA assay counts, usually `RNA`
- ATAC/chromatin assay counts, usually `ATAC` or `peaks`
- identical RNA and ATAC cell barcodes
- metadata columns:
  - `sample_id`
  - `cell_type`
  - optional `condition`

## Main workflow

```r
library(RegCompassR)

inputs <- rc_extract_seurat_v4(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)

pool_map <- rc_make_pools(
  inputs$meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition",
  target_size = 80,
  min_pool_size = 30,
  min_group_size = 30,
  seed = 1
)

layer1 <- rc_run_layer1_from_counts(
  gpr_table = gpr_table,
  rna_counts = inputs$rna_counts,
  pool_map = pool_map,
  stratum_col = "cell_type",
  tau = 0.20
)

C_raw <- layer1$C_raw
C_rel <- layer1$C_rel
reaction_confidence <- layer1$reaction_confidence
diagnostics <- layer1$minimal_diagnostics
```

## Mathematical defaults

### Pool expression

Layer 1 uses only:

\[
X^{RNA}_{g,p}=\log_2\left(1+\frac{count_{g,p}}{\sum_g count_{g,p}}\times 10^6\right)
\]

where \(count_{g,p}\) is the raw-count sum of gene \(g\) in pool \(p\).

Cell-level Pearson residuals, TF-IDF, scaled data, and imputed matrices are not valid main inputs because their averages are not equivalent to pseudobulk count normalization.

### Gene score

For each gene across pools:

\[
z_{g,p}=\frac{X_{g,p}-median(X_g)}{\max(MAD_\sigma, IQR/1.349, 0.05)}
\]

with clipping to \([-6,6]\), then:

\[
s_{g,p}=\sigma(z_{g,p})
\]

### GPR capacity

Promiscuous genes are downweighted by:

\[
s'_{g,p}=s_{g,p}/\sqrt{N_{rxn}(g)}
\]

For AND groups:

\[
w_{g,p}=\frac{\exp(-s'_{g,p}/\tau)}{\sum_h \exp(-s'_{h,p}/\tau)}
\]

\[
C_k(r,p)=\sum_g w_{g,p}s'_{g,p}
\]

Default \(\tau=0.20\). Smaller values such as 0.08 behave closer to a hard minimum; larger values move toward a mean-like AND. RegCompassR reports a hard-min sensitivity flag but keeps one main result.

For OR groups:

\[
C_{raw}(r,p)=\sum_k C_k(r,p)
\]

### Q95 calibration

Within each cell type:

\[
Q_r=\rho_n Q_{r,celltype}+(1-\rho_n)Q_{r,global}
\]

\[
\rho_n=\frac{n}{n+80}
\]

\[
C_{rel}(r,p)=\min\left(1,\frac{C_{raw}(r,p)}{Q_r+\epsilon}\right)
\]

`C_rel` is a relative reaction capacity potential, not enzyme activity and not a flux bound.

### Reaction confidence

When only RNA detection is available, reaction confidence defaults to:

\[
Conf_{r,p}=median(Detection_{g,p}:g\in GPR_r)\times(1-missing\_gpr\_gene\_fraction)
\]

If a user supplies a gene-level confidence matrix, the same missing-GPR penalty is applied.

## Minimal diagnostics

The main workflow returns:

- pool diagnostics:
  - `pool_id`
  - `sample_id`
  - `cell_type`
  - `condition`
  - `n_cells`
  - `low_power_pool`
  - `RNA_depth_mean`
  - `GPR_gene_detection_rate`
- Q95 diagnostics:
  - `q_stratum`
  - `q_global`
  - `rho_n`
  - `q95_low_power`
- GPR diagnostics:
  - `has_isoenzyme`
  - `has_multisubunit`
  - `missing_gpr_gene_fraction`
- tau sensitivity:
  - mean absolute difference from hard-min AND
  - `tau_sensitive_flag`
