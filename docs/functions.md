# Public functions in RegCompassR 1.8.1

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: prepare supported GEMs.
- `rc_make_medium_scenarios()`: create the shared extracellular medium.
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: run the complete GRN-first workflow.
- `rc_regcompass_step_grn()`: normalize single-cell RNA, compute cell-type-shared ATAC TF-IDF across conditions, and fit one Pando model per condition × cell type. Default `peak_cor = 0.01`.
- `rc_regcompass_step_metacells()`: build SuperCell2 metacells stratified only by condition. Default `gamma = 75`. Cell type is assigned afterwards from the dominant member-cell label. Exact dominant-label ties are rejected. The stage writes final root-level metadata, membership, composition and summary tables used downstream.
- `rc_regcompass_step_meta_modules()`: validate bidirectional GRN–metacell group coverage, write `grn_metacell_group_coverage.tsv.gz`, and convert condition × cell-type GRNs into complete-GPR core reactions and expanded biological meta-modules, followed by local FASTCORE completion.
- `rc_regcompass_step_layer1()`: integrate metacell RNA and ATAC evidence into reaction expression using the post hoc dominant cell-type label to select the matching GRN.
- `rc_regcompass_step_layer2()`: run directional COMPASS-like minimum-penalty scoring and maintain the persistent model cache.
- `rc_regcompass_step_results()`: assemble rankings and retain both condition-specific and global meta-module outputs.

Sample balancing APIs are not part of the workflow. `sample_col` is optional provenance only and does not affect cell selection, weighting, grouping, or graph construction. Cell type is not a metacell stratification variable.
