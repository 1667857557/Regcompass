# Select the RNA and ATAC reductions used for metacells

RegCompass Stage 2 constructs multimodal SuperCell metacells from two cell-level
embedding geometries:

- RNA reduction: `pca` dimensions `1:30` by default;
- ATAC reduction: `lsi` dimensions `2:30` by default;
- WNN neighbours: `k.knn = 30L` by default;
- graining target: `gamma = 30L` by default;
- base random seed: `12345L` by default.

For every broad cell type, all conditions are supplied jointly to one native
RNA+ATAC WNN graph. Condition is used only after Walktrap clustering to split
parent memberships into condition-pure final metacells. Fragment files, when
provided, are processed **after** this membership is fixed; they do not alter
the graph scope or introduce sample as a metacell stratum.

## Complete Stage 2 parameter example

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L,
    min_cells_per_stratum = 20L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L
  )
)
```

The implementation defaults are:

```r
rna_reduction = "pca"
rna_dims = 1:30
atac_reduction = "lsi"
atac_dims = 2:30
gamma = 30L
k.knn = 30L
kith = NULL
kernel = TRUE
seed = 12345L
min_cells_per_stratum = 20L
min_metacell_size = 1L
min_metacells_per_stratum = 1L
metacellNormalization = FALSE
avg.in.data = FALSE
verbose = FALSE
```

`min_cells_per_stratum` is checked on each final condition × cell-type coverage
stratum before graph construction. It does not cause separate condition graphs.
`min_metacell_size` only marks `low_power_metacell`; small final metacells are
not silently merged or removed.

The same fields can be supplied to `rc_run_regcompass()` or
`rc_run_regcompass_one_shot()` through `metacell_args`.

## Optional raw fragment input

Stage 2 supports two ATAC aggregation routes:

```text
fragment_files = NULL/FALSE
    -> aggregate the existing single-cell ATAC count matrix

fragment_files supplied
    -> build final SuperCell membership
    -> map single-cell fragment barcodes to final metacell IDs
    -> SuperCell::AggregateFragmentFile()
    -> optional MACS2 peak calling
    -> Signac::FeatureMatrix recount on metacell fragments
    -> metacell TF-IDF
```

For one fragment file, pass its path directly:

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  fragment_files = "/data/fragments.tsv.gz"
)
```

For multiple 10x fragment files when Seurat cell IDs are prefixed as
`sample_barcode`, use a named character vector whose names are the exact sample
prefixes:

```r
fragment_files <- c(
  Control_1 = "/data/Control_1/fragments.tsv.gz",
  Control_2 = "/data/Control_2/fragments.tsv.gz",
  Cre_1 = "/data/Cre_1/fragments.tsv.gz"
)
```

The prefix is used only to route cells to their source fragment file and recover
the raw barcode after the first `prefix_`. It is **not** used to stratify WNN or
metacell construction.

If cell IDs do not follow that convention, pass an explicit mapping data frame:

```r
fragment_map <- data.frame(
  fragment_file = c("/data/s1/fragments.tsv.gz", "/data/s2/fragments.tsv.gz"),
  object_cell = c("cell_name_in_Seurat_1", "cell_name_in_Seurat_2"),
  fragment_barcode = c("AAAC...-1", "AAAC...-1")
)
```

Each `object_cell` must belong to the exact Stage 2 cell set. Reused raw 10x
barcodes across different files are safe because mappings are validated within
each fragment file.

Fragment aggregation/recount controls are nested in
`metacell_args$fragment_args`:

```r
metacell_args = list(
  fragment_args = list(
    workers = NULL,
    rows_per_chunk = 10000000L,
    bgzip_path = NULL,
    tabix_path = NULL,
    process_n = 2000L,
    call_peaks = TRUE,
    macs2_path = NULL,
    effective_genome_size = NULL,
    peak_calling_args = list()
  )
)
```

When `call_peaks = TRUE`, MACS2 is run on the aggregated metacell fragment
files. With `call_peaks = FALSE`, the existing ATAC peak ranges are retained and
only counts are recomputed from fragments.

## Seed handling

`seed` is the reproducible base seed. SuperCell processes broad cell-type graph
groups in deterministic order and uses:

```text
seed_for_graph_group = seed + graph_group_index - 1
```

Repeated runs with the same cells, metadata, reductions, dimensions and seed
therefore use the same seed sequence. The Stage 2 input fingerprint records the
exact geometry and count inputs used for provenance.

## Use Harmony for the RNA geometry

A precomputed Harmony reduction can replace PCA without changing how membership
is constructed:

```r
stopifnot(
  "harmony" %in% names(A@reductions),
  "lsi" %in% names(A@reductions),
  ncol(SeuratObject::Embeddings(A[["harmony"]])) >= 30,
  ncol(SeuratObject::Embeddings(A[["lsi"]])) >= 30
)

step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  metacell_args = list(
    rna_reduction = "harmony",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L
  )
)
```

## What each geometry parameter changes

- `rna_reduction`: RNA embedding used by the multimodal WNN calculation.
- `rna_dims`: exact RNA coordinates supplied to SuperCell.
- `atac_reduction`: ATAC embedding used by the multimodal WNN calculation.
- `atac_dims`: exact ATAC coordinates supplied to SuperCell. LSI dimension 1 is
  excluded by default because it often tracks sequencing depth.
- `k.knn`: number of neighbours in the cell-type-scoped WNN graph.
- `gamma`: approximate cells-per-parent-metacell compression target.
- `seed`: deterministic base seed for graph-group construction.

Reduction names must exist in `A@reductions`, and the maximum requested
dimension must not exceed the number of columns in the corresponding embedding.

## What Harmony changes

Harmony changes the RNA neighbourhood geometry used by SuperCell to decide which
cells can be compressed together. It does not replace the RNA assay, alter RNA
counts or generate metacell expression values from Harmony coordinates. After
membership is determined, RegCompass always aggregates RNA counts. ATAC counts
come either from aggregation of the existing ATAC matrix or, when
`fragment_files` is supplied, from the fragment aggregation/recount route above.

Use Harmony only when it removes technical distortion without erasing the
biological condition contrast that RegCompass is intended to analyse. Never fit
separate Harmony coordinates by condition for this workflow: all conditions
must remain in one shared coordinate system within each broad cell type.

## Graph scope and condition purity

The graph boundary and purity boundary are different:

```text
graph boundary     = broad cell type
condition role     = post-clustering membership split
sample role        = fragment routing/provenance only when needed
```

Thus different broad cell types never share graph edges, while all conditions
within one broad cell type jointly determine the same local WNN geometry.

## Input fingerprint and restart behaviour

`step2$pooled$cache_contract` is retained as an input fingerprint and provenance
record. It contains:

- reduction names and selected dimensions;
- fingerprints of selected embeddings;
- ordered cells, broad cell types and conditions;
- RNA and ATAC count fingerprints used for WNN construction;
- grouped-WNN API and provenance fields;
- `gamma`, `k.knn`, `seed` and metacell thresholds.

When fragments are supplied, `step2$pooled$fragment_manifest` records the source
fragment files, generated metacell fragment files, indexes, metacell barcodes and
mapping route. The large fragment files themselves remain external files under
the Stage 2 output directory rather than being embedded in the RDS checkpoint.

Stage 2 no longer writes large sidecar RDS files for automatic cache recovery.
The authoritative restart point is `step_metacells.rds`. If Stage 2 inputs or
geometry change, run Stage 2 again and replace that checkpoint; there is no
`overwrite` control for a hidden sidecar cache.

## Verify the selected geometry and fragment route

```r
step2$pooled$cache_contract$analysis_args[c(
  "rna_reduction",
  "rna_dims",
  "atac_reduction",
  "atac_dims",
  "gamma",
  "k.knn",
  "seed",
  "min_cells_per_stratum",
  "min_metacell_size",
  "min_metacells_per_stratum"
)]

step2$pooled$input_design[c(
  "native_supercell_api",
  "graph_scope",
  "condition_scope",
  "membership_split_timing",
  "modality_weighting",
  "atac_aggregation_method",
  "fragment_files_supplied"
)]

step2$pooled$fragment_manifest
```

The embedding objects are fingerprints of the exact coordinates used to
construct the metacells.
