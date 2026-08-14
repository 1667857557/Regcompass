# Stage 2 metacell parameters

RegCompass Stage 2 uses two cell-level geometries:

- RNA: `pca`, dimensions `1:30` by default;
- ATAC: `lsi`, dimensions `2:30` by default.

For each broad cell type, all retained conditions are supplied to **one** multimodal WNN graph and **one** Walktrap hierarchy. Final condition-pure metacells are selected as condition-specific feasible cuts of that shared hierarchy. No condition-specific WNN graph and no post-hoc affinity repair are used.

## Common parameters

```r
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
```

`gamma` is a soft cells-per-metacell resolution target. `min_metacell_size` and `min_metacells_per_stratum` are hard constraints. Before graph construction RegCompass requires, for every condition × cell-type stratum,

```text
N_cells >= min_metacell_size * min_metacells_per_stratum
```

SuperCell then verifies whether the shared Walktrap hierarchy can satisfy those constraints. An infeasible hierarchy raises an error rather than emitting undersized metacells.

Use another RNA reduction only when it exists in the Seurat object and shares the same cell coordinate system, for example:

```r
metacell_args = list(rna_reduction = "harmony")
```

## Fragment recount

Supply `fragment_files` to `rc_regcompass_step_metacells()` only when ATAC counts should be rebuilt from fragments. Fragment routing occurs **after** final membership is fixed and does not change the WNN graph or hierarchy.

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = step1$params$requested_condition_col,
  celltype_col = step1$params$celltype_col,
  fragment_files = fragment_files,
  workers = workers
)
```

## Provenance

The canonical cell-to-metacell mapping is `step2$pooled$membership`. It retains `cell_id`, `metacell_id`, graph group, condition, shared hierarchy cut, community identifier, and partition policy. `step2$pooled$partition_diagnostics` records realized condition-specific cut statistics. The same final membership is used by RNA aggregation, ATAC aggregation/recount, and downstream Layer 1 projection.
