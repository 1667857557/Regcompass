# Metacell graph contract

RegCompass Stage 2 separates **graph scope** from **metacell purity**.

## Required architecture

For cell type \(t\), let

\[
I_t = \{i : c_i = t\}
\]

be all single cells of that cell type across every condition. RNA and ATAC
embeddings are standardized only within \(I_t\), using all conditions jointly:

\[
\widetilde Z^{(m)}_{t,ij} =
\frac{Z^{(m)}_{ij}-\mu^{(m)}_{t,j}}
     {s^{(m)}_{t,j}\sqrt{d_m}},
\qquad
\mu^{(m)}_{t,j}=\frac{1}{|I_t|}\sum_{i\in I_t}Z^{(m)}_{ij},
\]

where \(m\in\{RNA,ATAC\}\), \(d_m\) is the number of retained dimensions, and
zero-variance dimensions are set to zero. Dividing by \(\sqrt{d_m}\) gives each
modality block comparable total squared-distance weight.

The joint cell-type embedding is

\[
X_t = [\widetilde Z^{(RNA)}_t,\widetilde Z^{(ATAC)}_t].
\]

A separate kNN graph \(G_t=(I_t,E_t)\) is constructed for every cell type:

\[
E_t = \{(i,j): j\in kNN_{X_t}(i)\}.
\]

Therefore, for \(t\neq t'\), no edge can connect \(I_t\) and \(I_{t'}\). All
conditions inside \(I_t\) are present during distance calculation, neighbour
search, and graph clustering.

After graph clustering returns preliminary membership \(h_t(i)\), condition
purity is imposed by the interaction

\[
M(i)=\operatorname{interaction}(t,h_t(i),q_i),
\]

where \(q_i\) is condition. Condition does not define a graph and is not used to
standardize embeddings. It only splits a preliminary cluster when that cluster
contains multiple conditions.

This gives both required invariants:

1. **cell-type graph isolation**: \((i,j)\in E\Rightarrow c_i=c_j\);
2. **condition-pure metacells**: \(M(i)=M(j)\Rightarrow q_i=q_j\), while cells
   from different conditions may still be neighbours in the shared cell-type
   graph.

## SuperCell interface

RegCompass calls
`SuperCell::SCimplify_by_graph_group_from_embedding()` with:

- `cell.graph.group = cell_type`;
- `cell.split.condition = condition`;
- no sample-derived grouping field.

The formal provenance values are:

```text
native_supercell_api = SCimplify_by_graph_group_from_embedding
graph_scope = one_independent_graph_per_cell_type
condition_scope = all_conditions_joint_within_cell_type_graph
membership_split_timing = after_joint_graph_clustering
embedding_scaling = within_celltype_joint_condition_equal_modality_blocks
temporary_combined_stratum = FALSE
```

`cell.graph.group` determines which cells may enter the same graph.
`cell.split.condition` preserves condition-pure output membership after the
joint graph has been clustered.

## Why the previous interface was insufficient

Passing cell type as a post-clustering annotation to a single global graph can
prevent mixed metacells, but it does not prevent another cell type from changing
nearest-neighbour ranks, graph density, or clustering cuts before the split.
Conversely, splitting the input by condition before graph construction creates
condition-specific geometries and makes metacells less directly comparable.

The current contract avoids both errors: cell type is the graph boundary and
condition is a post-clustering purity boundary.

## Relationship to Pando modes

This Stage 2 contract is independent of automatic Stage 1 routing:

- `condition_grn` uses the canonical `pando_condition_grn_fit` output when at
  least two condition levels are present;
- `standard_pando` calls `Pando::infer_grn()` when condition metadata are absent
  or contain one level. No condition coefficients are calculated in this mode.

Both modes consume the same condition-pure metacell matrices and shared
cell-type graph provenance downstream.
