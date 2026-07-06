# RegCompassR

RegCompassR is a Layer 1 plus Layer 2 package for **sample-aware, multimodal-metacell-supported, GPR-aware reaction capacity potential** analysis from annotated Seurat/Signac RNA+ATAC objects. The main public workflow now starts from **SuperCell2.0 metacells** rather than random or embedding micropools.

## What the current code actually implements

```text
Annotated Seurat/Signac RNA+ATAC object
→ split by sample_id × condition × cell_type, plus optional state_col
→ SuperCell::SCimplify_for_Seurat() per eligible stratum
→ immediately save metacell_object.rds, membership.tsv.gz, metacell_metadata.tsv.gz,
  rna_counts.rds, atac_counts.rds, optional SuperCell fragment outputs, and QC files
→ import saved metacell counts/metadata
→ metacell RNA log2(CPM + 1)
→ robust gene z-score and sigmoid score through existing Layer 1 functions
→ GPR-aware reaction capacity, Q95 calibration, and reaction confidence
→ sample × condition × cell_type summaries from metacell-level scores
```

SuperCell2.0 is used for metacell construction and, when supported by the installed SuperCell version and fragment inputs, fragment aggregation. RegCompassR handles validation, count normalization, GPR-aware Layer 1 capacity/confidence, sample-level summaries, diagnostics, reporting, and selected-reaction Layer 2 utilities.

## Code-to-documentation alignment audit

| Tutorial claim | Current code status | Relevant function |
|---|---|---|
| Split by `sample_id × condition × cell_type` and optionally `state_col` | Implemented through `group_cols <- c(sample_col, condition_col, celltype_col, state_col)` after dropping missing grouping labels. | `rc_make_supercell2_metacells()` |
| Do not build metacells for tiny strata | Implemented: strata with fewer than `min_cells_per_stratum` only get a skipped QC file and no count matrices. | `rc_make_supercell2_metacells()` |
| Call SuperCell2.0 with RNA/ATAC assays, reductions, dims, gamma, optional labels, and optional fragment parameters | Implemented by constructing an argument list for `SuperCell::SCimplify_for_Seurat()`. | `rc_make_supercell2_metacells()` |
| Save metacell object, membership, metadata, RNA counts, ATAC counts, and QC immediately | Implemented for eligible strata. | `rc_make_supercell2_metacells()` |
| Save metacell fragment files | Passed through to SuperCell with `fragmentFiles`, `outputDirMcFragment`, `bgzip_path`, and `tabix_path` when `save_fragments = TRUE` and `fragment_files` is non-NULL. The exact fragment file names are produced by SuperCell; RegCompassR imports them with `fragments/*.tsv.gz`. | `rc_make_supercell2_metacells()`, `rc_import_supercell2_metacells()` |
| Run Layer 1 directly from raw metacell counts | Implemented. RNA counts are filtered for nonzero library size, converted to logCPM, and sent to the existing GPR/Q95/confidence engine. | `rc_run_layer1_from_metacells()` |
| Recompute peak-gene links on a metacell Signac object | Implemented as a wrapper around the existing Signac relinking helper. Requires Signac/qlcMatrix at runtime and a usable metacell Signac object. | `rc_recompute_metacell_peak_gene_links()` |
| Summarize biological samples, not treating metacells as replicates | Implemented as long-format medians over metacells within sample/condition/cell-type groups, with `n_metacells_used` and `single_metacell_group_flag`. | `rc_metacell_sample_summary()` |
| Old random/embedding pools are not the main public path | Implemented in NAMESPACE: old pool/pseudobulk APIs are not exported. The remaining `R/pseudobulk.R` helpers are compatibility internals used by legacy lower-level code and `.rc_aggregate_counts_by_membership()`. | `NAMESPACE`, `R/pooling.R`, `R/pseudobulk.R` |

## Expected input

RegCompassR expects an already annotated Seurat/Signac multiome object with:

- RNA raw UMI counts, usually assay `RNA`.
- ATAC raw peak counts, usually assay `ATAC` or `peaks`.
- Matching RNA and ATAC metacell/cell barcodes at the step where they are used.
- Metadata columns:
  - `sample_id`
  - `condition`
  - `cell_type`
  - optional `cell_state` or another state/subcluster column passed as `state_col`.
- Existing reductions needed by SuperCell2.0, commonly RNA PCA and ATAC LSI.
- Optional fragment files passed through `fragment_files` if metacell-level fragments are needed.

RegCompassR does **not** perform full QC, doublet removal, ambient correction, cell annotation, or WNN preprocessing. Complete those steps before running this package.

## Main workflow

```r
library(RegCompassR)
library(SuperCell)
library(Seurat)
library(Signac)

mc <- rc_make_supercell2_metacells(
  object = object,
  outdir = "RegCompassR_run/01_supercell2_metacells",
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
  fragment_files = fragment_files,
  save_metacell_object = TRUE,
  save_counts = TRUE,
  save_fragments = TRUE
)

# Optional: recompute metabolic peak-gene links on a saved metacell object.
# For multi-stratum runs, repeat or combine links across the desired metacell objects.
peak_gene_links <- rc_recompute_metacell_peak_gene_links(
  metacell_object = readRDS(mc$metacell_objects[[1]]),
  metabolic_genes = metabolic_genes,
  peak_assay = "ATAC",
  expression_assay = "RNA",
  out_file = "RegCompassR_run/02_peak_gene_links/metacell_peak_gene_links.tsv.gz"
)

layer1 <- rc_run_layer1_from_metacells(
  gpr_table = gpr_table,
  rna_metacell_counts = mc$rna_counts,
  metacell_meta = mc$metacell_meta,
  atac_metacell_counts = mc$atac_counts,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type",
  reaction_confidence_method = "gpr_aware",
  bootstrap = TRUE,
  B = 500
)

sample_summary <- rc_metacell_sample_summary(
  score_mat = layer1$C_rel,
  metacell_meta = layer1$metacell_meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)
```

If metacells have already been built and saved, skip construction and import them:

```r
mc <- rc_import_supercell2_metacells(
  metacell_dirs = Sys.glob("RegCompassR_run/01_supercell2_metacells/*"),
  rna_assay = "RNA",
  atac_assay = "ATAC",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  require_fragments = FALSE
)
```

## Saved metacell layout produced by the code

Eligible strata produce:

```text
RegCompassR_run/
  01_supercell2_metacells/
    sample_id=<sample>__condition=<condition>__cell_type=<cell_type>/
      metacell_object.rds              # only when save_metacell_object = TRUE
      membership.tsv.gz
      metacell_metadata.tsv.gz
      rna_counts.rds                   # only when save_counts = TRUE
      atac_counts.rds                  # only when save_counts = TRUE
      fragments/                       # SuperCell-created files, names may vary
        *.tsv.gz
        *.tsv.gz.tbi
      qc/
        metacell_qc.tsv.gz
        run_params.yaml                # only when yaml is installed
```

Strata below `min_cells_per_stratum` currently save only `qc/metacell_qc.tsv.gz` with `skipped = TRUE`; they are not imported into the returned count matrices.

`rc_import_supercell2_metacells()` returns:

```r
list(
  metacell_meta = data.frame(...),
  membership = data.frame(...),
  rna_counts = Matrix,
  atac_counts = Matrix or NULL,
  metacell_objects = character(),
  fragment_files = character(),
  diagnostics = data.frame(...)
)
```

## Public metacell API

- `rc_make_supercell2_metacells()` builds and saves SuperCell2.0 metacells per stratum.
- `rc_import_supercell2_metacells()` imports saved metacell directories.
- `rc_validate_metacell_inputs()` validates count/metadata alignment and RNA/ATAC metacell column agreement.
- `rc_build_metacell_metadata()` builds one-row-per-metacell metadata from a membership table.
- `rc_filter_empty_metacells()` removes zero-library metacells before logCPM.
- `rc_metacell_detection()` converts raw metacell counts to gene-by-metacell binary detection.
- `rc_atac_metacell_logcpm()` filters ATAC peaks by metacell detection and computes logCPM.
- `rc_run_layer1_from_metacells()` runs RegCompass Layer 1 from raw metacell RNA counts and optional ATAC/link evidence.
- `rc_recompute_metacell_peak_gene_links()` recomputes metabolic peak-gene links on a metacell Signac object.
- `rc_metacell_sample_summary()` summarizes metacell reaction scores at sample/condition/cell-type level.
- `rc_metacell_diagnostics()` and `rc_write_metacell_report()` generate basic metacell QC outputs.

## Human-GEM GPR tables and metabolic peak-gene links

```r
hg <- rc_download_humangem_gpr_table(
  destdir = "data/Human-GEM",
  ref = "main",
  gene_format = "symbol"
)

gpr_table <- hg$gpr_table
metabolic_genes <- hg$metabolic_genes
```

For multiome confidence, prefer peak-gene links recomputed on metacell objects when metacell fragments are available. If no fragments are available, users may supply an external `peak_gene_links` table with `peak_id`, `gene`, and `weight` columns.

## Guardrails and limitations

- Do not construct formal-analysis metacells across samples, conditions, or major cell types.
- Do not use integrated/corrected expression as RegCompass input; Layer 1 expects raw metacell count sums followed by logCPM.
- Do not attach original single-cell fragment files directly to a metacell object. Use SuperCell-generated metacell fragments when fragment-level Signac operations are needed.
- Metacells are denoised analysis units, not biological replicates. Formal comparisons should use biological samples.
- `gamma = 75`, RNA dims `1:30`, ATAC dims `2:30`, `min_cells_per_stratum = 100`, and `min_metacell_size = 20` are defaults, not universal biological constants. Run sensitivity analyses for key conclusions.
- `cell_state` should be used as `state_col` or `label_col` only when each sample/state stratum has enough cells.
- Layer 2 should be run on Layer 1-selected reactions to avoid unnecessary large LPs.
- In this implementation, `rc_run_layer1_from_metacells()` still calls existing lower-level Layer 1 helpers whose argument names use `pool_*` internally; this is a compatibility detail, not the public workflow terminology.

## Parallel execution

Metacell construction/import loops, reaction-wise GPR capacity, Q95 bootstrap diagnostics, and model-level Layer 2 tasks use `rc_parallel_lapply()` where implemented. Control workers with:

```r
options(RegCompassR.workers = 4)
# or before starting R:
# Sys.setenv(REGCOMPASS_WORKERS = "4")
```

Set `options(RegCompassR.workers = 1)` or pass `BPPARAM = FALSE` for deterministic sequential execution.
