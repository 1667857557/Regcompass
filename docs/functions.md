# Public function reference

Complete argument definitions are in the Rd help pages.

## Workflow

### `rc_run_regcompass_one_shot()`

Species-aware complete workflow. Main arguments: `object`, `outdir`, `genome`, `species`, optional `gem`, `gem_version`, `gem_source`, `pfm`, `fragment_files`, `medium_scenario`, `medium_scenarios`, `workers`, `progress`, and `...` forwarded to `rc_run_regcompass()`.

### `rc_run_regcompass()`

Complete workflow with an explicit GEM. Main arguments: `object`, `gem`, `outdir`, `genome`, `species`, `condition_col`, `celltype_col`, `cell_type`, `rna_assay`, `atac_assay`, `fragment_files`, `pando_args`, `target_rsq_threshold`, `metacell_args`, `meta_module_args`, `layer1_args`, `medium_scenarios`, `model_mode`, `layer2_args`, `workers`, `progress`.

## Restartable stages

### `rc_regcompass_step_grn()`

Stage 1 Pando GRN. `pando_args$min_cells = 500L` by default. Routing is automatic per broad cell type: multi-condition groups use condition common-dictionary ridge; one-condition groups use standard Pando. `target_rsq_threshold = 0.05` is the RegCompass final-fit target R² quality gate.

### `rc_regcompass_step_metacells()`

Stage 2 shared-WNN metacells. Key defaults: RNA `pca` 1:30, ATAC `lsi` 2:30, `gamma = 30L`, `k.knn = 30L`, `min_cells_per_stratum = 20L`, `min_metacell_size = 1L`, `min_metacells_per_stratum = 1L`. If `min_metacell_size > 1L`, supply `min_merge_affinity`. Optional raw ATAC fragments use `fragment_files`.

### `rc_regcompass_step_meta_modules()`

Stage 3 reaction catalogue and meta-modules. Optional customization uses `meta_module_args`.

### `rc_regcompass_step_layer1()`

Stage 4 reaction evidence. Main controls: `gpr_and_method`, `gene_half_saturation`, `workers`, `progress`.

### `rc_regcompass_step_layer2()`

Stage 5 structural model and directional scoring. `model_mode = "meta_module_gem"` uses CORDA2 by default. `layer2_args` accepts current Rd-documented structural/scoring controls.

### `rc_regcompass_step_results()`

Assemble final reaction results.

### `rc_regcompass_step_target_union()`

Rescore selected database-linked reaction targets in existing Stage 5 structural models.

## GEM and medium

### `rc_prepare_gem()`

Prepare the supported human or mouse GEM.

### `rc_make_medium_scenarios()`

Create predefined or custom medium tables. Presets: `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, `custom`. Built-in challenge concentrations are metadata unless an explicit flux/bound assumption is supplied. See [medium-presets.md](medium-presets.md).

## Parallel configuration

### `rc_parallel_config()`

Inspect the resolved worker/backend configuration. All stages share the top-level `workers` cap.

## Post-analysis

- `rc_test_condition_reactions()` — test selected reaction targets between conditions.
- `rc_plot_condition_reaction()` — plot one reaction target.
- `plot_top_celltype_reaction_rank()` — plot top directional reactions within a cell type.
- `rc_build_reaction_annotations()` / `rc_attach_reaction_annotations()` — build/attach reaction annotations.
- `rc_select_gene_reactions()` / `rc_plot_condition_gene_reactions()` — select and inspect GPR-associated reactions for genes.
