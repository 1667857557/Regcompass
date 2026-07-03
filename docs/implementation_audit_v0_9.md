# RegCompassR v0.9 implementation audit

This audit summarizes the current code-level scope after simplifying RegCompassR to a Layer 1 reaction-capacity workflow. It supersedes the older v0.6-v0.7 selected-demand QP planning notes.

## Current package scope

RegCompassR now focuses on a reproducible, diagnostic Layer 1 workflow:

- validate annotated Seurat v4 RNA+ATAC objects and extract raw counts with `rc_validate_seurat()` / `rc_extract_inputs()` plus the compatibility aliases `rc_validate_seurat_v4()` / `rc_extract_seurat_v4()`;
- create sample-aware, optional condition/state-aware micropools with `rc_make_pools()` and `rc_make_pool_seed_replicates()`;
- aggregate raw RNA counts by pool with `rc_pseudobulk_counts()`, remove empty pools with `rc_filter_empty_pools()`, and normalize to pool-level `log2(CPM + 1)` with `rc_logcpm()`;
- compute robust gene scores, GPR-aware reaction capacities, Q95 calibration, reaction confidence, and sensitivity diagnostics through `rc_run_layer1_from_counts()` / `rc_run_layer1_capacity()`;
- download Human-GEM and prepare a RegCompass-compatible GPR table plus the corresponding metabolic GPR gene set;
- optionally compute ATAC-supported multiome gene confidence from pooled ATAC accessibility and curated or externally regenerated peak-gene links;
- summarize, export, and report pool-level results at sample-aware levels.

## Removed QP planning layer

The following modules are intentionally out of scope in the current simplified package and their source files now contain removal placeholders:

- baseline GEM/QP construction;
- selected-demand QP sweeps;
- reaction-selection planning for QP workloads;
- regulator-ranking/causal-driver layers.

Tutorials should therefore describe `C_rel` as relative reaction capacity potential and should not present it as a hard flux bound or as direct flux inference.

## Input and metadata checks

The current input API requires a Seurat object with paired RNA and ATAC assays containing the same cell barcodes. Metadata checks cover required sample and cell-type columns plus optional condition, batch, and state columns. `rc_check_metadata()` and `rc_write_input_summary()` provide human-readable input summaries including sample/cell-type counts, missing metadata counts, optional condition-by-batch tables, and state-source records.

## Pooling and pseudobulk checks

Pooling remains sample-aware: cells are never pooled across samples, and optional condition/cell-type/state grouping columns define independent pooling strata. `rc_drop_na_grouping()` removes cells with missing grouping values before pool construction. `rc_check_pseudobulk_mapping()` provides spot checks that pseudobulk columns are consistent with pool membership.

## Layer 1 capacity and diagnostics

The core capacity calculation keeps raw counts until after pool aggregation. Gene scores use robust z-scores over pool-level logCPM values followed by a sigmoid transformation. GPR AND rules use the default Boltzmann minimum-biased average with `tau = 0.20`; OR rules sum isoenzyme-group capacities. Q95 calibration uses continuous shrinkage toward global Q95 values and can report bootstrap uncertainty.

The workflow returns capacity matrices and diagnostics including Q95 power classes, GPR gene coverage, hard-min/tau/promiscuity/AND-method sensitivity summaries, long-form capacity tables, parsed GPR rules, pool metadata, and the source of reaction confidence.


## Human-GEM and metabolic peak-gene links

`rc_download_humangem_gpr_table()` downloads a Human-GEM GitHub archive, reads `model/genes.tsv`, `model/reactions.tsv`, and `model/Human-GEM.yml`, converts Human-GEM gene identifiers to symbols by default, and returns a RegCompass-compatible `gpr_table`, `metabolic_genes`, raw reaction rules, and source annotation tables.

`rc_metabolic_gpr_genes()` extracts the metabolic GPR gene set from any parsed or tabular GPR input. That gene set is intended as the target for regenerating metabolic peak-gene links in an external multiome workflow such as Signac `LinkPeaks(genes.use = metabolic_genes)`. The package does not run peak-gene link inference internally; `rc_run_layer1_from_counts()` accepts a supplied `peak_gene_links` table, filters it to the GPR metabolic genes, and then computes multiome confidence from pooled ATAC accessibility and link weights.

## Multiome confidence

When pooled ATAC counts and curated peak-gene links are supplied, the wrapper filters links to genes present in the GPR set, computes pooled ATAC logCPM, derives RNA and ATAC percentiles within the selected stratum, estimates null-corrected RNA/ATAC concordance, applies Fisher shrinkage to positive association evidence, and combines nonnegative components into gene confidence. Single-pool strata produce undefined percentiles rather than artificially high confidence.

## Sample-aware summaries and reports

`rc_sample_aggregate()` aggregates pool-level matrices to biological sample × annotated cell-type medians. `rc_sample_summary()` returns long-form sample/cell-type/condition summaries with median and IQR. Export helpers write sample matrices and long pool-level tables, and `rc_write_report_md()` generates a compact Markdown diagnostic report.
