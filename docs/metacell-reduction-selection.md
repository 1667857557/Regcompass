# Stage 2 metacell geometry and repair

RegCompass Stage 2 constructs multimodal SuperCell metacells from cell-level RNA and ATAC embeddings. The default geometry is:

- RNA: `pca`, dimensions `1:30`;
- ATAC: `lsi`, dimensions `2:30`;
- WNN neighbours: `k.knn = 30L`;
- parent graining target: `gamma = 30L`;
- base seed: `12345L`.

The canonical membership sequence is fixed:

```text
broad cell type
  -> all retained conditions jointly
  -> one native SuperCell RNA+ATAC WNN
  -> Walktrap parent clustering
  -> split parent memberships by condition
  -> optional small-metacell repair on that exact original shared WNN
  -> final condition-pure membership
  -> RNA/ATAC aggregation and all downstream RegCompass calculations
```

Condition is never used to build a separate WNN. Small-metacell repair never rebuilds a condition-specific graph.

## Routine Stage 2 call

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  workers = 10L
)
```

Default metacell controls are:

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
min_merge_affinity = NULL
unresolved_small_policy = "error"
metacellNormalization = FALSE
avg.in.data = FALSE
verbose = FALSE
```

With the default `min_metacell_size = 1L`, no repair is required and the condition-split membership is retained. `min_cells_per_stratum` is checked for each condition × broad-cell-type stratum before graph construction; it does not create separate condition graphs.

## Enforce a minimum final metacell size

When `min_metacell_size > 1L`, an explicit `min_merge_affinity` must also be supplied:

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  metacell_args = list(
    min_metacell_size = 10L,
    min_metacells_per_stratum = 3L,
    min_merge_affinity = 0.05,
    unresolved_small_policy = "error"
  ),
  workers = 10L
)
```

The numeric threshold above is only an example; RegCompass does not supply a hidden biological default for `min_merge_affinity`.

For provisional metacells `M` and `N`, the repair affinity is computed from the exact original cell-type shared WNN after non-negative symmetrization:

\[
W \leftarrow \frac{W + W^T}{2},
\]

\[
A(M,N)=\frac{\sum_{i\in M,j\in N}W_{ij}}
{\sqrt{\operatorname{vol}(M)\operatorname{vol}(N)}},\qquad
\operatorname{vol}(M)=\sum_{i\in M}d_i.
\]

A small metacell can merge only with a metacell from the same condition and the same broad cell type, and only when `A(M,N) >= min_merge_affinity`. Repair is deterministic: smaller groups are considered first, stable IDs break ordering ties, and affinities are recomputed after each merge. A merge is never allowed to reduce a stratum below `min_metacells_per_stratum`.

The necessary cell-count feasibility condition is

\[
N_{stratum} \ge min\_metacell\_size\times min\_metacells\_per\_stratum.
\]

This is not sufficient: the original WNN topology or the affinity threshold can still leave an undersized metacell. With `unresolved_small_policy = "error"`, Stage 2 stops. With `"keep"`, the unresolved group is retained and explicitly reported.

The final `membership_table` contains the repaired `metacell_id` and the original `parent_metacell_id`. Parent IDs are provenance only: one repaired final metacell may contain cells from multiple original Walktrap parents.

## Downstream membership contract

The final repaired membership is the single canonical membership used by every downstream calculation:

- SuperCell RNA/ATAC aggregation;
- fragment aggregation/recount when raw ATAC fragments are supplied;
- quantitative RNA `mean(single-cell CPM)` aggregation;
- Pando regulatory projection `beta * mean(TF) * mean(ATAC)`;
- Layer 1 and Layer 2 reaction evidence.

RegCompass does not retain an alternative pre-repair membership for scoring. The original WNN is returned transiently by `SuperCell::SCimplify_for_Seurat(return.graph = TRUE)` only inside the grouped builder, is used for repair, and is discarded immediately afterward.

## Use another RNA reduction

A precomputed Harmony reduction can replace PCA without changing the membership algorithm:

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  metacell_args = list(
    rna_reduction = "harmony",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30
  ),
  workers = 10L
)
```

Reduction names must exist in `A@reductions`, and each requested dimension must exist in the corresponding embedding. LSI dimension 1 is excluded by default because it often tracks sequencing depth. Harmony changes only RNA neighbourhood geometry; it does not replace the RNA assay or determine metacell expression values. Do not fit separate condition-specific embeddings for this workflow: all conditions within a broad cell type must remain in a common coordinate system.

## Optional raw ATAC fragment input

Without fragments, Stage 2 aggregates the existing single-cell ATAC count matrix. When fragments are supplied, RegCompass first fixes the final repaired membership, then maps fragment barcodes to those final IDs, aggregates fragments, optionally calls peaks, and recounts the metacell peak matrix.

For one file:

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  outdir = "RegCompass_steps/02_condition_metacells",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  fragment_files = "/data/fragments.tsv.gz",
  workers = 10L
)
```

For prefixed Seurat cell IDs, a named vector can route cells to source files:

```r
fragment_files <- c(
  Control_1 = "/data/Control_1/fragments.tsv.gz",
  Control_2 = "/data/Control_2/fragments.tsv.gz",
  Cre_1 = "/data/Cre_1/fragments.tsv.gz"
)
```

The prefix is fragment provenance only and is not a WNN/metacell stratum. For other naming schemes, pass an explicit `fragment_map` with `fragment_file`, `object_cell`, and `fragment_barcode` columns.

Fragment controls are nested under `metacell_args$fragment_args`; worker count remains the top-level `workers` argument:

```r
metacell_args = list(
  fragment_args = list(
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

## Provenance and restart contract

`step2$pooled$cache_contract` fingerprints the inputs and algorithmic controls that define membership, including:

- selected RNA/ATAC reductions and dimensions;
- selected embedding fingerprints;
- ordered cells, conditions, and broad cell types;
- RNA/ATAC count fingerprints;
- `gamma`, `k.knn`, `seed`;
- `min_metacell_size`, `min_metacells_per_stratum`, `min_merge_affinity`;
- `unresolved_small_policy`;
- repair algorithm/version and WNN symmetrization rule.

Repair diagnostics are available in `step2$pooled$repair_diagnostics`; unresolved groups, if permitted, are in `step2$pooled$unresolved_small_metacells`. `parent_metacell_id` in the membership table preserves the original Walktrap parent.

The authoritative restart artifact is `step_metacells.rds`. If membership-defining inputs or controls change, Stage 2 must be regenerated.
