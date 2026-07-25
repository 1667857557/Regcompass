# Public functions in RegCompassR 1.8.4

## Setup and complete runs

- `rc_prepare_gem()`, `rc_prepare_human2_gem()`, `rc_prepare_mouse_gem()`: load a bundled pinned model or explicitly download/rebuild a release.
- `rc_bundled_gem_manifest()`: inspect installed Human-GEM/Mouse-GEM release, checksum, citation, and license metadata.
- `rc_download_species_gem()`: lower-level official repository download/parse path.
- `rc_parallel_config()`: inspect OS-specific backend resolution.
- `rc_make_medium_scenarios()`: create one shared medium table; see [medium presets](medium-presets.md).
- `rc_run_regcompass()` and `rc_run_regcompass_one_shot()`: execute the complete significant-Pando-target workflow with progress, timing, automatic backend selection, `upstream_workers = 6L`, and `layer2_workers = 30L`.

The complete workflow exposes two worker counts. `upstream_workers` covers Pando inference and Layer 1. `layer2_workers` covers medium-specific union-GEM construction, global FASTCORE completion, and directional LP scoring. Stage 3 performs biological catalogue construction without FASTCORE.

## Inspectable stages

- `rc_regcompass_step_grn()`: fit condition-by-cell-type Pando models for all Human-GEM GPR genes present in the RNA assay. Human analyses default to the union of Pando's hg38 phastCons and SCREEN ccRE data objects unless `pando_initiate_args$regions` is supplied.
- `rc_regcompass_step_metacells()`: condition-level, cell-type-guided SuperCell2 metacells. RNA PCA is the default geometry; an existing Harmony reduction can be selected through `metacell_args$rna_reduction` and `metacell_args$rna_dims`.
- `rc_regcompass_step_meta_modules()`: summarize significantly supported metabolic target genes, map complete-GPR cores, perform subsystem/KEGG/Reactome/master-Rhea expansion, and deduplicate reaction IDs into a merged catalogue. No target projection, connected-component analysis, FASTCORE, or GEM construction occurs here.
- `rc_regcompass_step_layer1()`: integrated RNA+ATAC reaction support.
- `rc_regcompass_step_layer2()`: first construct one medium-specific union GEM with global FASTCORE, then cache it and run directional LP scoring; or use shared full-GEM scoring when `model_mode = "full_gem"`.
- `rc_regcompass_step_results()`: rankings, reaction annotations, evidence provenance, and condition contrasts.
- `rc_regcompass_step_target_union()`: map selected original core reactions through shared KEGG, Reactome, or master-Rhea identifiers and score mapped non-core reactions in the exact final Stage 5 union GEMs.

## Stage 3 evidence and catalogue fields

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

`meta_module_args` accepts only `subsystem_table`, `expansion_mode`, and `max_iterations`. `layer1_args` accepts only `regulatory_alpha`, `tau`, and `gene_half_saturation`.

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
- `rc_test_condition_reactions()`: descriptive same-target comparisons across conditions within cell type.
- `rc_select_gene_reactions()`: select scored reactions by GPR gene.
- `rc_plot_condition_reaction()` and `rc_plot_condition_gene_reactions()`: annotated condition plots.

Sample balancing is not part of the canonical workflow. Metacell-level comparisons are descriptive pseudo-observation analyses rather than automatic biological-replicate inference.

## Tutorials

- [Portable execution, bundled GEMs, progress, timing, and worker cleanup](portable-execution.md)
- [Metacell PCA/Harmony reduction selection](metacell-reduction-selection.md)
- [Level 1: quick start](tutorial-01-quick-start.md)
- [Level 2: true stepwise workflow](tutorial-02-stepwise-audit.md)
- [Level 3: restart and diagnostics](tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](tutorial-05-condition-differential-analysis.md)
