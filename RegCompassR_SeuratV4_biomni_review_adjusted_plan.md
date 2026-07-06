# RegCompassR SuperCell2.0 Metacell Tutorial and Implementation Plan

This document supersedes the earlier micropool/pseudobulk tutorial. It is aligned with the current code in `R/metacell.R`, `NAMESPACE`, and `README.md`.

## 1. Current design target

```text
Annotated Seurat/Signac RNA+ATAC object
→ sample_id × condition × cell_type split, plus optional state_col
→ SuperCell2.0 metacell construction per eligible stratum
→ immediate persistence of metacell object, membership, metadata, raw RNA/ATAC counts, optional fragments, and QC
→ RegCompass Layer 1 from metacell raw RNA counts, with optional ATAC peak-gene confidence
→ sample × condition × cell_type reaction summaries
→ selected-reaction Layer 2 utilities
```

The public workflow is metacell-first. Random/embedding micropool construction is not exported as a main API.

## 2. Inputs expected by the current implementation

Required for `rc_make_supercell2_metacells()`:

```text
object: Seurat object
RNA assay: default RNA
ATAC assay: default ATAC
metadata columns:
  sample_id
  condition
  cell_type
optional metadata:
  state_col
  label_col
reductions:
  RNA reduction name, default pca
  ATAC reduction name, default lsi
```

Optional fragment support:

```text
fragment_files
bgzip_path
tabix_path
```

The code passes these fragment arguments to `SuperCell::SCimplify_for_Seurat()`. Fragment output names are determined by SuperCell and are imported by globbing `fragments/*.tsv.gz`.

## 3. Metacell construction function

Implemented public function:

```r
rc_make_supercell2_metacells(
  object,
  outdir,
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  state_col = NULL,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  rna_reduction = "pca",
  atac_reduction = "lsi",
  rna_dims = 1:30,
  atac_dims = 2:30,
  gamma = 75,
  min_cells_per_stratum = 100,
  min_metacell_size = 20,
  label_col = NULL,
  fragment_files = NULL,
  bgzip_path = "bgzip",
  tabix_path = "tabix",
  save_metacell_object = TRUE,
  save_counts = TRUE,
  save_fragments = TRUE,
  overwrite = FALSE,
  BPPARAM = NULL
)
```

Behavior implemented in code:

1. Require package `SuperCell` at runtime.
2. Validate Seurat class and required metadata columns.
3. Drop cells with missing split labels using `rc_drop_na_grouping()`.
4. Split cells by `sample_col`, `condition_col`, `celltype_col`, and optional `state_col`.
5. Write a skipped QC record for strata below `min_cells_per_stratum`.
6. For eligible strata, subset the Seurat object and call `SuperCell::SCimplify_for_Seurat()`.
7. Pass fragment arguments only when `save_fragments = TRUE` and `fragment_files` is non-NULL.
8. Rename metacell columns with a stratum prefix and `MC###` suffix.
9. Save `metacell_object.rds`, `membership.tsv.gz`, `metacell_metadata.tsv.gz`, `rna_counts.rds`, `atac_counts.rds`, `qc/metacell_qc.tsv.gz`, and `qc/run_params.yaml` when enabled and available.
10. Import the saved eligible strata with `rc_import_supercell2_metacells()` and return merged counts/metadata.

## 4. Saved outputs

Eligible strata produce this code-backed layout:

```text
<outdir>/
  sample_id=<sample>__condition=<condition>__cell_type=<cell_type>/
    metacell_object.rds
    membership.tsv.gz
    metacell_metadata.tsv.gz
    rna_counts.rds
    atac_counts.rds
    fragments/
      *.tsv.gz
      *.tsv.gz.tbi
    qc/
      metacell_qc.tsv.gz
      run_params.yaml
```

Important correction from the old tutorial: skipped strata currently save only a QC file and are not included in returned matrices.

## 5. Importing existing metacells

Implemented public function:

```r
rc_import_supercell2_metacells(
  metacell_dirs,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  require_fragments = FALSE
)
```

The importer reads saved metadata/membership/counts, validates RNA/ATAC column agreement when ATAC counts exist, checks fragment indexes if `require_fragments = TRUE`, and returns:

```r
list(
  metacell_meta,
  membership,
  rna_counts,
  atac_counts,
  metacell_objects,
  fragment_files,
  diagnostics
)
```

## 6. Layer 1 from metacells

Implemented public function:

```r
rc_run_layer1_from_metacells(
  gpr_table,
  rna_metacell_counts,
  metacell_meta,
  atac_metacell_counts = NULL,
  peak_gene_links = NULL,
  metacell_seurat = NULL,
  recompute_peak_gene_links = FALSE,
  metabolic_genes = NULL,
  stratum_col = "cell_type",
  promiscuity_mode = "sqrt",
  and_method = "boltzmann",
  tau = 0.20,
  reaction_confidence_method = "gpr_aware",
  bootstrap = TRUE,
  B = 500,
  BPPARAM = NULL
)
```

Implemented flow:

```text
raw RNA metacell counts
→ remove zero-library metacells
→ rc_logcpm()
→ rc_metacell_detection()
→ optional ATAC metacell logCPM + peak-gene confidence
→ rc_run_layer1_capacity()
→ return C_raw/C_rel/confidence plus metacell metadata and metacell logCPM/detection
```

The lower-level Layer 1 engine still has historical `pool_*` argument names internally; the public metacell wrapper converts `metacell_id` to an internal `pool_id` column only for compatibility.

## 7. Peak-gene links on metacells

Implemented public function:

```r
rc_recompute_metacell_peak_gene_links(
  metacell_object,
  metabolic_genes = NULL,
  peak_assay = "ATAC",
  expression_assay = "RNA",
  genome = NULL,
  distance = 5e5,
  min_cells = 10,
  out_file = NULL,
  gpr_table = NULL,
  ...
)
```

This wraps the existing Signac helper and writes a gzipped link table if `out_file` is supplied. It requires Signac and qlcMatrix at runtime. Provide either `metabolic_genes` or `gpr_table` so the metabolic gene set can be determined.

## 8. Sample-level summary

Implemented public function:

```r
rc_metacell_sample_summary(
  score_mat,
  metacell_meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = NULL
)
```

This returns a long table with per-reaction medians across metacells and diagnostic columns:

```text
group_id
reaction_id
value
n_metacells_used
n_cells_used
low_power_group_flag
single_metacell_group_flag
```

## 9. Diagnostics and report

Implemented public functions:

```r
rc_metacell_diagnostics()
rc_write_metacell_report()
```

Current diagnostics include RNA depth, RNA detected genes, ATAC detected peaks, and GPR gene detection rate when the required matrices/genes are supplied.

## 10. What was removed from the tutorial

The old tutorial sections on `rc_make_pools()`, pool seed replicates, random pooling, embedding pooling, and user-facing pseudobulk pool APIs are no longer valid as a main workflow and have been removed from this plan. The package still contains some legacy lower-level pseudobulk helpers for compatibility, but they are not exported as the public metacell workflow.
