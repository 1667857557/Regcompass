# Function reference

This page lists the supported public API. Complete argument definitions are in the corresponding Rd help pages. Equations and quantitative definitions are maintained only in [mathematical-model.md](mathematical-model.md).

## Complete workflow

### `rc_run_regcompass_one_shot()`

Prepare species defaults when needed and run the complete workflow.

Main arguments: `object`, `outdir`, `genome`, `species`, optional `gem`, `gem_version`, `gem_source`, `pfm`, `fragment_files`, `medium_scenario`, `medium_scenarios`, `workers`, `progress`, plus `...` forwarded to `rc_run_regcompass()`.

### `rc_run_regcompass()`

Run the complete workflow with an explicit GEM and stage argument bundles.

Main arguments: `object`, `gem`, `outdir`, `genome`, `species`, `condition_col`, `celltype_col`, `cell_type`, `rna_assay`, `atac_assay`, `fragment_files`, `pando_args`, `metacell_args`, `meta_module_args`, `layer1_args`, `medium_scenarios`, `model_mode`, `layer2_args`, `workers`, and `progress`.

## Restartable stages

### `rc_regcompass_step_grn()`

Infer Stage 1 regulatory evidence. `pando_args$min_cells` defaults to `500L`; inference options are supplied through `pando_args$pando_infer_args`. Multi-condition cell types use Pando multi-task ridge. Standard Pando now defaults to `method = "ridge"`, which is the K=1 specialization of the same ridge solver; `method = "glm"` remains an explicit compatibility option. Ridge routes keep one cell type resident at a time and use target-level workers under the global worker cap. Each Pando target worker receives only the target-specific RNA, ATAC, motif and edge-dictionary slices required for that target, and worker/batch temporaries are released after completion.

### `rc_regcompass_step_metacells()`

Construct Stage 2 multimodal metacells. The public `metacell_args$min_cells_per_stratum` default is `500L`. Optional raw fragments are supplied through `fragment_files`.

### `rc_regcompass_step_meta_modules()`

Build the Stage 3 reaction catalogue and cell-type reaction meta-modules. Optional customization is supplied through `meta_module_args`.

### `rc_regcompass_step_layer1()`

Project RNA and regulatory evidence to reactions. Public controls are `gpr_and_method`, `gene_half_saturation`, `workers`, and `progress`.

### `rc_regcompass_step_layer2()`

Build the Stage 5 structural model and score directional targets. `layer2_args` accepts `model_params`, `omega`, `target_direction`, `solver`, and `flux_threshold`. `model_mode = "meta_module_gem"` uses CORDA2 by default; supplementary structural routes are selected explicitly through the current Rd-documented controls.

### `rc_regcompass_step_results()`

Assemble reaction annotations, evidence, rankings, metacell-level comparisons, and available condition contrasts.

### `rc_regcompass_step_target_union()`

Select directly database-linked targets and rescore them in existing Stage 5 structural models without rebuilding the structural model.

## GEM and medium

### `rc_prepare_gem()`

Prepare the supported species GEM. Important arguments are `species`, `version`, `source`, `cache_dir`, `save_rds`, `force_download`, and `allow_latest`. Human defaults to Human-GEM `2.0.0`; mouse defaults to Mouse-GEM `1.8.0`.

Species-specific preparation wrappers and low-level download/asset helpers are not public workflow APIs. Use `rc_prepare_gem()` for both species.

### `rc_make_medium_scenarios()`

Create built-in or custom medium tables. Supported presets are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. See [medium-presets.md](medium-presets.md).

## Parallel configuration

### `rc_parallel_config()`

Inspect the platform-resolved parallel configuration without starting workers. The workflow itself uses one top-level `workers` cap. Ridge GRNs do not nest cell-type and target pools: cell types are processed sequentially and the available worker budget is reused for target-level Pando work. Target payloads are trimmed before dispatch and are released promptly after each completed worker batch.

## Post analysis

### `rc_test_condition_reactions()`

Compare fixed reaction/direction/medium targets between conditions within cell type.

### `rc_plot_condition_reaction()`

Plot one selected reaction target across conditions.

### `plot_top_celltype_reaction_rank()`

Plot top directional reaction targets within one cell type.

### `rc_build_reaction_annotations()` / `rc_attach_reaction_annotations()`

Build reaction annotations and attach them to a result.

### `rc_select_gene_reactions()` / `rc_plot_condition_gene_reactions()`

Select GPR-associated reactions for specified genes and optionally test/plot condition differences.
