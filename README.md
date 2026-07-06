# RegCompassR

RegCompassR is a Layer 1 plus Layer 2 tool for **sample-aware multimodal metacell-supported, GPR-aware reaction capacity potential** analysis from annotated Seurat/Signac RNA+ATAC objects. The main analysis unit is now a **SuperCell2.0 metacell**, not a random or embedding micropool.

The package intentionally keeps the main analysis narrow:

```text
Annotated Seurat/Signac RNA+ATAC object
→ sample × condition × cell_type SuperCell2.0 multimodal metacells
→ save metacell object, membership, raw RNA/ATAC metacell counts, and metacell fragments
→ metacell-level log2(CPM + 1)
→ robust z-score by gene across metacells
→ sigmoid gene score
→ sqrt promiscuity correction
→ GPR capacity:
   AND = Boltzmann minimum-biased average, tau = 0.20
   OR  = sum across isoenzyme groups
→ cell-type Q95 continuous shrinkage
→ GPR-aware multiome/RNA reaction confidence
→ sample × cell_type reaction summary
→ selected-subnetwork Layer 2 penalty LP on Layer 1-selected reactions
```

SuperCell2.0 is responsible for sample-aware RNA+ATAC metacell construction and fragment aggregation. RegCompassR is responsible for GPR-aware reaction capacity, multiome confidence, Q95 calibration, sample-level summaries, reporting, and Layer 2 selected-subnetwork modeling.

## Expected input

RegCompassR expects an already annotated Seurat/Signac multiome object with:

- RNA raw UMI counts, usually assay `RNA`.
- ATAC raw peak counts in a ChromatinAssay, usually assay `ATAC` or `peaks`.
- Identical RNA and ATAC cell barcodes.
- Fragment files attached to the ATAC assay or supplied as `fragment_files` when building metacells.
- Metadata columns:
  - `sample_id`
  - `condition`
  - `cell_type`
  - optional `cell_state` / subcluster labels.
- Existing RNA PCA and ATAC LSI reductions for SuperCell2.0.

RegCompassR does **not** perform full QC, doublet removal, ambient correction, cell annotation, or WNN preprocessing. Those steps should be completed before calling RegCompassR.

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
  rna_assay = "RNA",
  atac_assay = "ATAC",
  rna_reduction = "pca",
  atac_reduction = "lsi",
  rna_dims = 1:30,
  atac_dims = 2:30,
  gamma = 75,
  min_cells_per_stratum = 100,
  min_metacell_size = 20,
  fragment_files = fragment_files,
  save_metacell_object = TRUE,
  save_counts = TRUE,
  save_fragments = TRUE
)

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
  atac_metacell_counts = mc$atac_counts,
  metacell_meta = mc$metacell_meta,
  peak_gene_links = peak_gene_links,
  stratum_col = "cell_type",
  reaction_confidence_method = "gpr_aware",
  bootstrap = TRUE,
  B = 500
)

sample_summary <- rc_metacell_sample_summary(
  score_mat = layer1$C_rel,
  metacell_meta = mc$metacell_meta,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  condition_col = "condition"
)
```

## Saved metacell layout

`rc_make_supercell2_metacells()` saves every stratum immediately so interrupted runs do not lose expensive metacell or fragment work:

```text
RegCompassR_run/
  01_supercell2_metacells/
    sample_id=<sample>__condition=<condition>__cell_type=<cell_type>/
      metacell_object.rds
      membership.tsv.gz
      metacell_metadata.tsv.gz
      rna_counts.rds
      atac_counts.rds
      fragments/
        metacell_fragments.tsv.gz
        metacell_fragments.tsv.gz.tbi
      qc/
        metacell_qc.tsv.gz
        run_params.yaml
  02_peak_gene_links/
  03_regcompass_layer1/
  04_sample_summary/
  05_reports/
```

The required saved artifacts are:

1. `membership.tsv.gz`: original cell ID to metacell ID mapping.
2. `metacell_object.rds`: SuperCell2.0 Seurat/Signac metacell object.
3. `rna_counts.rds` and `atac_counts.rds`: raw metacell counts used by RegCompass.
4. `metacell_metadata.tsv.gz`: sample, condition, cell type, metacell size, and low-power flags.
5. `fragments/*.tsv.gz` and `.tbi`: metacell-level fragment files when fragments are supplied.
6. `qc/run_params.yaml`: gamma, dims, assays, fragment paths, and software settings when `yaml` is installed.

Existing outputs can be loaded with `rc_import_supercell2_metacells()`.

## API overview

- `rc_make_supercell2_metacells()` splits by `sample_id × condition × cell_type` (and optional state), calls SuperCell2.0, saves outputs, and returns merged metacell counts and metadata.
- `rc_import_supercell2_metacells()` reads saved metacell directories and validates count/metadata alignment.
- `rc_validate_metacell_inputs()` checks metacell raw counts, ATAC/RNA column agreement, and metadata uniqueness.
- `rc_filter_empty_metacells()`, `rc_metacell_detection()`, and `rc_atac_metacell_logcpm()` provide metacell-level normalization helpers.
- `rc_recompute_metacell_peak_gene_links()` recomputes metabolic Signac peak-gene links on metacell objects.
- `rc_run_layer1_from_metacells()` runs RegCompass Layer 1 directly from raw metacell RNA/ATAC counts.
- `rc_metacell_sample_summary()` summarizes metacell reaction scores at the biological `sample × condition × cell_type` level.
- `rc_metacell_diagnostics()` and `rc_write_metacell_report()` create metacell QC summaries.

The old random/embedding micropool construction functions are no longer exported as a main workflow. Historical pseudobulk helpers remain internal fallbacks for membership-based aggregation and compatibility with lower-level code, but new analyses should start from SuperCell2.0 metacell raw counts.

## Human-GEM GPR tables and metabolic peak-gene links

RegCompassR can download Human-GEM GPR rules and convert them into the simple `reaction_id`, `and_group_id`, `gene` table used by Layer 1:

```r
hg <- rc_download_humangem_gpr_table(
  destdir = "data/Human-GEM",
  ref = "main",
  gene_format = "symbol"
)

gpr_table <- hg$gpr_table
metabolic_genes <- hg$metabolic_genes
```

For multiome confidence, recompute peak-gene links on the metacell object rather than on original single cells whenever metacell fragments are available.

## Biological and engineering guardrails

- Do not construct formal-analysis metacells across samples, conditions, or major cell types.
- Do not use integrated/corrected expression as RegCompass input; Layer 1 uses raw metacell count sums followed by logCPM.
- Do not attach original single-cell fragment files directly to a metacell object; use SuperCell2.0 metacell fragments.
- Treat metacells as denoised analysis units, not biological replicates. Formal comparisons remain sample-level.
- Use `gamma = 75`, RNA dims `1:30`, ATAC dims `2:30`, `min_cells_per_stratum = 100`, and `min_metacell_size = 20` as defaults, then run gamma sensitivity (`50`, `100`) for key results.
- Use `cell_state` mainly as a semi-supervised label or sensitivity analysis unless every sample × state has enough cells.
- Run Layer 2 only for Layer 1-selected reactions to keep LP runtimes tractable.

## Parallel execution

RegCompassR parallelizes reaction-wise GPR capacity, Q95 bootstrap diagnostics, metacell construction/import loops, and model-level Layer 2 work when `BiocParallel` is installed. Control workers with:

```r
options(RegCompassR.workers = 4)
# or before starting R:
# Sys.setenv(REGCOMPASS_WORKERS = "4")
```

Set `options(RegCompassR.workers = 1)` or pass `BPPARAM = FALSE` for deterministic sequential execution.
