# Metacell graph contract

RegCompass Stage 2 separates **graph scope** from **metacell purity**.

## Required architecture

For broad cell type \(t\), define

\[
I_t = \{i : celltype_i=t\}.
\]

All cells in \(I_t\), across every condition, are supplied together to one
multimodal RNA+ATAC WNN construction:

\[
G_t = \operatorname{WNN}\left(
Z^{RNA}_{I_t}, Z^{ATAC}_{I_t}
\right).
\]

The native SuperCell/Seurat WNN implementation learns adaptive RNA and ATAC
modality weights within that cell type. Different broad cell types never share a
single-cell graph:

\[
(i,j)\in E_t \Rightarrow celltype_i=celltype_j=t.
\]

Walktrap clustering on \(G_t\) produces parent membership \(h_t(i)\). Condition
is not used to define graphs, choose neighbours, or fit separate condition
geometries. It is applied only after clustering:

\[
M(i)=\operatorname{interaction}\left(t,h_t(i),condition_i\right).
\]

Therefore final metacells satisfy both invariants:

1. **cell-type graph isolation:** graph edges never connect different broad cell types;
2. **condition purity:** cells in one final metacell have one condition, while all conditions participated jointly in the parent WNN graph.

Small condition-split metacells are retained. Stage 2 marks them with
`low_power_metacell` according to `min_metacell_size`; it does not merge or
remove them.

## SuperCell interface

RegCompass calls:

```r
SuperCell::SCimplify_by_graph_group(
  seurat = object,
  cell.graph.group = object[[celltype_col]][, 1],
  cell.split.condition = object[[condition_col]][, 1],
  assay = c(rna_assay, atac_assay),
  reduction = list(rna_reduction, atac_reduction),
  dims = list(rna_dims, atac_dims),
  k.knn = k.knn,
  gamma = gamma
)
```

Canonical provenance values are:

```text
native_supercell_api = SCimplify_by_graph_group
graph_group_argument = cell.graph.group
condition_argument = cell.split.condition
graph_method = multimodal_WNN
graph_scope = one_independent_WNN_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_WNN_graph_clustering
modality_weighting = adaptive_WNN_within_cell_type
temporary_combined_stratum = FALSE
sample_metadata = not_used_or_retained
```

RNA PCA and ATAC LSI reductions must already exist in a shared coordinate system
for all conditions. They must not be recomputed separately by condition.

## Count aggregation

The grouped WNN call determines membership. RegCompass then passes the exact
condition-pure membership to `SCimplify_for_Seurat()` in membership mode to
aggregate RNA and ATAC counts. It does not rebuild the WNN graph during count
aggregation.

## Relationship to Pando

Stage 1 Pando and Stage 2 SuperCell use the same broad cell-type grouping but
serve different roles:

- Pando jointly fits condition-aware TF–peak–gene effects within each cell type;
- SuperCell jointly builds the WNN geometry within each cell type and splits
  membership by condition after clustering;
- RegCompass aggregates cell-level OOF regulatory projections using the exact
  final membership.

No `sample` column is required by the canonical workflow.
