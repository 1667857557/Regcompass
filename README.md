# RegCompassR

RegCompassR is an R package scaffold for RegCompass-Multiome: a multiome-supported, GPR-aware, sample-aware framework for reaction capacity potential and selected network-constrained feasibility analysis.

## v0.1 scope

The v0.1 implementation focuses on the input layer for already annotated Seurat v4/Signac single-cell multiome objects. It does not rerun clustering, WNN, or metabolic modeling.

Implemented functions:

- `rc_validate_seurat()` checks that a Seurat object contains the requested RNA assay, ATAC assay, required sample and cell-type metadata, optional condition/batch metadata, and optional embedding.
- `rc_extract_inputs()` validates the object and extracts RNA assay data, ATAC assay data, cell metadata, and an optional embedding into a plain R list.

## Expected input

A Seurat v4 object should contain:

- RNA assay: commonly `RNA` or `SCT`
- ATAC/chromatin assay: commonly `ATAC` or `peaks`
- metadata columns such as `sample_id`, `cell_type`, optional `condition`, and optional `batch`
- optional reductions such as WNN, UMAP, PCA, LSI, or Harmony

## Example

```r
library(RegCompassR)

rc_validate_seurat(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)

inputs <- rc_extract_inputs(
  object,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  celltype_col = "cell_type"
)
```
