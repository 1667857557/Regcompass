# Public functions in RegCompassR 1.9.1

## Setup and complete runs

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: load a bundled pinned model or explicitly download/rebuild a release.
- `rc_bundled_gem_manifest()`: inspect installed Human-GEM/Mouse-GEM release, checksum, citation, and license metadata.
- `rc_download_species_gem()`: lower-level official repository download/parse path.
- `rc_parallel_config()`: inspect OS-specific backend resolution.
- `rc_make_medium_scenarios()`: create one shared medium table; see [medium presets](medium-presets.md).
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: execute the complete significant-Pando-target workflow with progress, timing, automatic backend selection, `upstream_workers = 6L`, and `layer2_workers = 30L`.

Public runner arguments are arranged by processing sequence: shared inputs, Stage 1 Pando, Stage 2 metacells, Stage 3 meta-modules, Stage 4 Layer 1, Stage 5 Layer 2, and execution controls.

The complete workflow exposes two worker counts. `upstream_workers` covers Pando inference and Layer 1. `layer2_workers` covers medium-specific union-GEM construction, global FASTCORE completion, and directional LP scoring. Stage 3 performs biological catalogue construction without FASTCORE.

## Inspectable stages

- `rc_regcompass_step_grn()`: fit one shared-design Pando model per cell type with independently estimated, directly comparable condition coefficients for all GEM GPR genes present in the RNA assay. When `pfm` is omitted, Pando's bundled `motifs` data object is used. Human defaults to phastCons plus SCREEN ccRE regions; mouse defaults to phastCons regions only. Explicit `pando_initiate_args$regions` overrides either default.
- `rc_regcompass_step_metacells()`: condition-level, cell-type-guided SuperCell2 metacells. RNA PCA dimensions 1:30, ATAC LSI dimensions 2:30, and `seed = 12345L` are the defaults; an existing Harmony reduction can replace PCA through `metacell_args$rna_reduction` and `metacell_args$rna_dims`.
- `rc_regcompass_step_meta_modules()`: summarize significantly supported metabolic target genes, map complete-GPR cores, perform one fixed ordered subsystem/KEGG/Reactome/master-Rhea expansion pass, and deduplicate reaction IDs into a merged catalogue. No target projection, connected-component analysis, FASTCORE, or GEM construction occurs here.
- `rc_regcompass_step_layer1()`: integrated RNA+ATAC reaction support with COMPASS-compatible GPR-AND aggregation.
- `rc_regcompass_step_layer2()`: first construct one medium-specific union GEM with global FASTCORE, then cache it and run directional LP scoring; or use shared full-GEM scoring when `model_mode = "full_gem"`.
- `rc_regcompass_step_results()`: rankings, reaction annotations, evidence provenance, and condition contrasts.
- `rc_regcompass_step_target_union()`: map selected original core reactions through shared KEGG, Reactome, or master-Rhea identifiers and score mapped non-core reactions in the exact final Stage 5 union GEMs.

## Stage 1 evidence controls

The canonical Pando fit and evidence defaults are:

```r
pando_args = list(
  min_cells = 20L,
  min_abs_estimate = 0,
  min_model_rsq = 0.1,
  pando_infer_args = list(
    method = "shared_design_independent",
    candidate_screen = "condition_union",
    condition_mix = 1,
    condition_weight = "equal",
    scale = TRUE
  )
)
```

Active condition coefficients define the supported-gene table. Explicit
condition-minus-reference coefficients and the stored Pando transform feed the
Layer 1 regulatory projection.

## Stage 2 geometry and reproducibility

```r
metacell_args = list(
  rna_reduction = "pca",
  rna_dims = 1:30,
  atac_reduction = "lsi",
  atac_dims = 2:30,
  gamma = 30L,
  seed = 12345L,
  min_cells_per_stratum = 100L,
  min_metacell_size = 20L,
  min_metacells_per_stratum = 2L
)
```

For ordered condition strata, the SuperCell2 seed is `seed + stratum_index - 1`. These fields and the selected embedding fingerprints participate in cache validation.

## Stage 3 evidence and catalogue fields

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$meta_module_summary$expansion_policy
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

The recorded expansion policy is:

```text
single_ordered_annotation_pass
```

`meta_module_args` accepts only an optional `subsystem_table`. Expansion is always one ordered pass:

```text
core subsystem
→ KEGG/Reactome reaction equivalence
→ master-Rhea reaction equivalence
```

`layer1_args` accepts only `regulatory_alpha`, `gpr_and_method`, and `gene_half_saturation`. `gpr_and_method` accepts `"min"`, `"median"`, or `"mean"` and defaults to `"min"`.

The merged object is a biological reaction catalogue, not a GEM. It contains no medium constraints or FASTCORE support.

## Stage 5 union-GEM configuration

```r
layer2_args$model_params <- list(
  completion_time_limit = 600,
  fastcore_epsilon = 1e-4,
  max_support_reactions = 2000,
  strict = TRUE
)
```

`completion_time_limit` applies exclusively to FASTCC/FASTCORE union-GEM construction. The directional scoring API has no `time_limit` parameter and begins only after the structural model has been completed and cached.

Each row of `step5$model_cache_summary` identifies one final medium-specific union GEM and records its file checksum, reaction counts, FASTCORE support count, build strategy, and completion stage.

The optional second-pass scoring function uses these exact cached files. It validates the checksum and medium identity, does not rebuild a GEM, does not rerun FASTCORE, and has no scoring timeout control.

Stages validate workflow parameters, GEM fingerprints, stage classes, and ordered metacell IDs before accepting an upstream object. Current condition-metacell artifacts also contain a cache contract covering ordered cells, labels, assays, selected PCA/Harmony and LSI embeddings, construction labels, and analysis parameters. Every public stage returns a timing table and writes `step_timing.tsv`.

## Interpretation and plotting

- `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()`: reaction names, formulas, GPRs, and evidence classes.
- `rc_test_condition_reactions()`: descriptive same-target, same-direction comparisons across conditions within cell type.
- `rc_report_condition_directions()`: preserve forward/reverse LP results, diagnose numerically indistinguishable directions, and add non-additive best-direction and directional-asymmetry summaries without claiming net flux.
- `rc_select_gene_reactions()`: select scored reactions by GPR gene.
- `rc_plot_condition_reaction()` and `rc_plot_condition_gene_reactions()`: annotated condition plots.

Sample balancing is not part of the canonical workflow. Metacell-level comparisons are descriptive pseudo-observation analyses rather than automatic biological-replicate inference.

## Tutorials

- [Portable execution, bundled GEMs, progress, timing, and worker cleanup](portable-execution.md)
- [Metacell PCA/Harmony reduction selection](metacell-reduction-selection.md)
- [Direction-aware final reporting](direction-aware-condition-reporting.md)
- [Level 1: quick start](tutorial-01-quick-start.md)
- [Level 2: true stepwise workflow](tutorial-02-stepwise-audit.md)
- [Level 3: restart and diagnostics](tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](tutorial-05-condition-differential-analysis.md)
