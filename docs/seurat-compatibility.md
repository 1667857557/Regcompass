# Seurat v4/v5 compatibility

## Supported runtime profiles

| Profile | SeuratObject | Seurat | Signac | Default | Input assay policy |
|---|---:|---:|---:|---|---|
| `seurat_v4_default` | `>= 4.1.4, < 5` | `>= 4.4.0, < 5` | `>= 1.11.0, < 1.14.0` | yes | v3 `Assay` plus Signac 1.x `ChromatinAssay` |
| `seurat_v5_compatible` | `>= 5, < 6` | `>= 5, < 6` | `>= 1.12.0, < 2` | no | v3 `Assay` or joinable `Assay5`, plus Signac 1.x `ChromatinAssay` |

Seurat and SeuratObject must use the same major version. Signac 2.x introduces a separate `ChromatinAssay5` object model and is intentionally outside the current compatibility contract.

The package `DESCRIPTION` retains the exact default pins:

```text
SeuratObject 4.1.4
Seurat 4.4.0
Signac 1.11.0
```

These pins define the canonical reproducibility environment; compatibility with later versions does not change the default data model.

## Recommended Seurat v5 object policy

When Seurat v5 is required for other parts of an analysis, create v3-style assays for RegCompass inputs whenever possible:

```r
options(Seurat.object.assay.version = "v3")
```

This keeps the same single `counts` and `data` matrices used by the default Seurat v4 workflow while allowing the surrounding session to use Seurat v5.

## Existing Assay5 objects

Seurat v5 can store one logical assay as several layers, commonly:

```text
counts.sample1
counts.sample2
data.sample1
data.sample2
```

RegCompass requires one aligned feature-by-cell matrix for each semantic layer. The canonical Stage 1 and Stage 2 entry points therefore operate on a local working copy and apply the following policy:

1. detect the assay class rather than inferring storage from the installed package version;
2. retain a single exact `counts` or `data` layer unchanged;
3. join multiple `counts.*` or `data.*` layers with `SeuratObject::JoinLayers()`;
4. reject an ambiguous assay containing both an exact layer and split variants;
5. validate named feature and cell dimensions after joining;
6. record the installed stack, object version, assay classes, storage type, and joined layer names.

The original object supplied by the caller is not rewritten. R copy-on-modify semantics are used, and the prepared object is passed only through the current stage.

## Matrix access policy

All RegCompass assay access is routed through one internal compatibility layer:

| Operation | Seurat v4 / v3 Assay | Seurat v5 / Assay5 |
|---|---|---|
| read counts | `GetAssayData(slot = "counts")` fallback | `LayerData(layer = "counts")` |
| read normalized data | `GetAssayData(slot = "data")` fallback | `LayerData(layer = "data")` |
| write normalized data | `SetAssayData(slot = "data")` fallback | `LayerData(..., layer = "data") <- value` |
| split-layer resolution | not applicable | `JoinLayers()` on the working copy |

The implementation does not inspect the installed SeuratObject version to guess assay storage. A v3 `Assay` remains a v3 `Assay` even when it is used inside a Seurat v5 session.

## Provenance

Prepared objects contain:

```r
object@misc$regcompass_seurat_compatibility
```

The record includes:

```text
default_input_profile
installed_stack_profile
package_versions
object_version
assay_classes
assay_storage
joined_layers
```

Stage 2 also copies this record into:

```r
step2$params$seurat_compatibility
```

## Unsupported or rejected cases

RegCompass stops rather than choosing silently when:

- Seurat and SeuratObject major versions differ;
- a required assay or `counts` matrix is absent;
- an Assay5 contains both `counts` and `counts.*`, or both `data` and `data.*`;
- split layers cannot be joined;
- feature or cell identifiers are missing, duplicated, or inconsistent;
- Signac 2.x `ChromatinAssay5` is installed.

## Validation matrix

The default CI job remains pinned to the exact Seurat v4 stack. A separate Seurat v5 job verifies:

- package loading under the v5 profile;
- v3-style `Assay` access under Seurat v5;
- Assay5 split-layer detection and joining;
- normalized-data read/write behavior;
- preservation of feature and cell order;
- compatibility provenance.

Algorithmic workflow tests continue to run in the default Seurat v4 environment so that adding v5 support does not silently change canonical results.
