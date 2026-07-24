# Public functions in RegCompassR 1.8.1

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: prepare supported GEMs.
- `rc_make_medium_scenarios()`: create the shared extracellular medium. `physiologic` resolves to species-specific plasma; culture formulations, nutrient sensitivity scenarios, technical baselines, and custom media are described in [Predefined extracellular medium scenarios](medium-presets.md).
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: run the complete GRN-first workflow. On Linux, `upstream_workers` controls Stage 1 Pando groups, Stage 3 local FASTCORE completion, and Stage 4 GPR capacity; `layer2_workers` controls Stage 5 LP tasks. Use `parallel_backend = "multicore"` for explicit forked execution.
- `rc_regcompass_step_grn()`: normalize single-cell RNA, align assay data by cell ID, compute cell-type-shared ATAC TF-IDF across conditions, and fit one Pando model per condition × cell type. Default `peak_cor = 0.01`. Use `parallel = TRUE` with `BiocParallel::MulticoreParam` on Linux, while keeping Pando's inner `parallel = FALSE`.
- `rc_regcompass_step_metacells()`: build SuperCell2 metacells with condition as the only hard stratum and automatically pass `celltype_col` to SuperCell2 as the pre-aggregation label to reduce annotated cell-type mixing. Default `gamma = 75`. Cell type and purity are audited afterwards from member-cell labels. The canonical stage has no separate label-column argument and is not controlled by the workflow `BPPARAM`.
- `rc_regcompass_step_meta_modules()`: validate bidirectional GRN–metacell group coverage, convert GRNs into complete-GPR core reactions and expanded biological modules, then run local FASTCORE completion. Parallel completion is configured through `layer1_args$local_fastcore_args`, using `parallel`, `workers`, `backend`, or `BPPARAM`.
- `rc_regcompass_step_layer1()`: integrate metacell RNA and ATAC evidence into reaction expression. GPR/reaction-capacity work accepts `parallel` and `BPPARAM`.
- `rc_regcompass_step_layer2()`: preflight the selected LP solver, construct the persistent structural-model cache, and distribute shared-model × metacell LP tasks through `parallel` and `BPPARAM`. The default HiGHS solver is a required dependency.
- `rc_regcompass_step_results()`: assemble rankings and retain both condition-specific and global meta-module outputs. Stage 6 now also creates `reaction_catalog` and `reaction_evidence`, including names, formulas, GPR genes, and RNA-versus-multiome provenance.
- `rc_build_reaction_annotations()`: reconstruct formal reaction formulas from GEM stoichiometry and summarize GPR genes and evidence provenance.
- `rc_attach_reaction_annotations()`: add the reaction catalog and evidence tables to results generated before this annotation layer was available.
- `rc_test_condition_reactions()`: compare one shared-GEM reaction target across two or more conditions within each selected cell type. It reports Kruskal-Wallis omnibus tests, pairwise Wilcoxon tests, BH-adjusted P values, median shifts, rank-biserial/common-language effects, Cohen's d, reaction names/formulas, GPR genes, and evidence classes.
- `rc_select_gene_reactions()`: select scored reactions by one or more metabolic genes while retaining the full Boolean GPR rule and condition-by-cell-type evidence provenance.
- `rc_plot_condition_reaction()`: draw one selected `reaction × direction × medium × cell type` across conditions as a boxplot with one point per metacell, an omnibus subtitle, pairwise significance brackets, and a biological caption containing the tested formula, genes, and evidence source. Requires the suggested package `ggplot2`.
- `rc_plot_condition_gene_reactions()`: select significant reaction-direction targets by metabolic genes, rank them once using the complete multiple-testing family, and return a named collection of annotated boxplots.
- `rc_available_workers()` and `rc_default_bpparam()`: detect workers and construct an automatic, serial, socket, or multicore backend.

Sample balancing APIs are not part of the workflow. `sample_col` is optional provenance only and does not affect cell selection, weighting, grouping, or graph construction. Cell type is not a metacell stratification variable.

Metacell-level condition tests describe within-dataset separation and are not automatically biological-replicate-level treatment inference. See [Condition-associated reaction statistics](condition-reaction-statistics.md) for the score definition, reaction annotations, gene-driven selection, multiple-testing scope, plotting examples, and inference boundary.

## Tutorial level by API surface

- [Level 1](tutorial-01-quick-start.md) uses `rc_run_regcompass_one_shot()` with explicit Linux worker counts.
- [Level 2](tutorial-02-stepwise-audit.md) uses all six `rc_regcompass_step_*()` functions with `MulticoreParam` examples and input/output audit gates.
- [Level 3](tutorial-03-advanced-restart.md) covers worker allocation, restart, serial debugging, solver, medium, and model-scope diagnostics.
- [Condition statistics, reaction annotations, and plots](condition-reaction-statistics.md) covers condition tests, formula/GPR annotation, RNA-versus-multiome evidence, gene-based reaction selection, and plot collections.
