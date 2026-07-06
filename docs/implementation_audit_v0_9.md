# RegCompassR Implementation Audit: Metacell-first Code Alignment

This audit reflects the current metacell-first implementation rather than the historical pool/micropool design.

## Public workflow now implemented

- `rc_make_supercell2_metacells()` builds SuperCell2.0 metacells per `sample_id × condition × cell_type` stratum, with optional `state_col` and `label_col`.
- `rc_import_supercell2_metacells()` imports saved metacell outputs and validates count/metadata alignment.
- `rc_validate_metacell_inputs()`, `rc_build_metacell_metadata()`, `rc_filter_empty_metacells()`, `rc_metacell_detection()`, and `rc_atac_metacell_logcpm()` provide metacell-level input and normalization utilities.
- `rc_run_layer1_from_metacells()` runs Layer 1 from raw metacell RNA counts and optional ATAC/peak-gene evidence.
- `rc_recompute_metacell_peak_gene_links()` recomputes metabolic peak-gene links on metacell Signac objects when Signac dependencies and fragments are available.
- `rc_metacell_sample_summary()`, `rc_metacell_diagnostics()`, and `rc_write_metacell_report()` summarize and report metacell-level results.

## Important code-backed limitations

- Strata below `min_cells_per_stratum` are skipped and currently save QC only.
- Fragment files are delegated to SuperCell. RegCompassR passes fragment-related arguments and imports files matching `fragments/*.tsv.gz`; it does not guarantee fixed fragment file names.
- The Layer 1 internals still use some historical `pool_*` argument names. The public wrapper hides this by creating an internal `pool_id` alias from `metacell_id`.
- Legacy pseudobulk functions remain in `R/pseudobulk.R` for compatibility and internal fallback use, but are not exported in `NAMESPACE`.
- Automated R tests could not be executed in the current environment because `Rscript` is unavailable.

## NAMESPACE audit

The public namespace exports the metacell APIs and no longer exports `rc_make_pools()`, `rc_make_pool_seed_replicates()`, `rc_pseudobulk_counts()`, `rc_filter_empty_pools()`, `rc_build_pool_metadata()`, `rc_pool_detection()`, `rc_atac_pool_logcpm()`, or `rc_pool_diagnostics()`.

## Test audit

The repository now includes `tests/testthat/test_metacell.R`, which checks:

1. metacell metadata construction and validation;
2. Layer 1 execution from toy raw metacell counts;
3. sample-summary output fields for metacell diagnostics.

Historical pool/pseudobulk/diagnostics tests were removed because those APIs are no longer the documented main path.
