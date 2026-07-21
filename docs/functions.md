# Public functions in RegCompassR 1.8.0

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: prepare supported GEMs.
- `rc_make_medium_scenarios()`: create the shared extracellular medium.
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: run the complete GRN-first workflow.
- `rc_regcompass_step_grn()`: normalize single-cell RNA, compute cell-type-shared ATAC TF-IDF across conditions, and fit one Pando model per condition × cell type. Default `peak_cor = 0.01`.
- `rc_regcompass_step_metacells()`: build condition × cell type SuperCell2 metacells. Default `gamma = 75`.
- `rc_regcompass_step_meta_modules()`: convert condition-specific GRNs into complete-GPR core reactions and expanded biological meta-modules, followed by local FASTCORE completion.
- `rc_regcompass_step_layer1()`: integrate metacell RNA and ATAC evidence into reaction expression.
- `rc_regcompass_step_layer2()`: run directional COMPASS-like minimum-penalty scoring.
- `rc_regcompass_step_results()`: assemble rankings and the final result.

Sample balancing APIs are not part of the workflow. `sample_col` is optional provenance only and does not affect cell selection, weighting, grouping, or graph construction.
