- Union meta-modules, union GEMs, FASTCORE completion, model caches, and directional vmax reuse are now partitioned by cell type; only conditions within the same cell type share a structural model.
# RegCompassR 2.2.7

- Requires Pando 1.6.3 native ABI 6 and its explicit high-dimensional memory
  contract. The master process and up to two supplied BiocParallel workers now
  execute the registered native numerical/allocation self-test before Stage 1.
- Routes condition-aware Stage 1 through compact diagnostics, bounded target
  batches, and resumable checkpoints under the actual Stage 1 output directory.
  User memory budgets and other explicit engine controls remain authoritative.
- Records dense versus sparse matrix-free target counts, largest predictor and
  nonzero sizes, PCG diagnostics, and estimated peak bytes without changing the
  ConditionGRNFit coefficient, mask, transform, nested-CV, or OOF contracts.
- Resolves automatic downstream media only from recorded human or mouse GEM
  provenance. GEMs without species provenance now require an explicit
  biologically justified medium instead of silently selecting a retired
  technical boundary mode.

# RegCompassR 2.2.1

- Added structured stage progress for Stage 1. Console events now report phase, percent, elapsed time, cell types, conditions, metabolic targets, outer/inner folds, lambda path size, solver route, refit/validation and OOF status. Every run writes `step_progress.tsv`, including when console messages are disabled.
- Stage 1 now enables Pando target-level verbose output when RegCompass progress is enabled, while retaining the existing public API and fail-fast behavior.
- The RegCompass Pando route now skips unused exact motif-hit positions and
  retains only the binary peak-by-motif incidence matrix. Interrupted and
  failed stages immediately report their last phase and point to the progress
  and timing logs for diagnosis.
- Standard Pando now has the same detailed Stage 1 observability as the
  condition-aware route, including per-cell-type candidate, motif, fit,
  extraction and artifact events. Errors and interrupts print their original
  condition message together with the active phase before unwinding.

- Restored the documented biological medium scenarios `normal_human_plasma`,
  `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`,
  `low_lactate`, and `low_glutamine`, with separate composition and challenge
  provenance.
- Replaced RPMI/DMEM and lower-tier serum surveys as canonical composition
  sources. Human composition now uses HPLM from *Cell* 2017 and updated HPLM
  from *Cell Metabolism* 2021; Plasmax from *Science Advances* 2019 is retained
  as independent validation and is not numerically averaged with HPLM.
- Anchored `mouse_plasma` to absolute mouse plasma and interstitial-fluid
  metabolomics in *Nature* 2026. Unsupported mouse components are omitted
  rather than inherited from human HPLM; glucose, lactate, and glutamine retain
  limited quantitative secondary values.
- All five culture challenges now use the identical authoritative HPLM 2017/2021
  background and override only the named glucose, lactate, or glutamine
  concentration. Outputs mark the construction as
  `authoritative_HPLM_background_plus_named_nutrient_override`.
- Corrected `low_glutamine` to the Methods-defined 0.5 mM condition from Visagie
  et al. rather than the prior 0.05 mM value.
- Retained reaction-level and metabolite-level user-defined media through
  `scenario = "custom"` or `scenario = NULL`, including mixed built-in and
  custom scenario runs.
- Kept technical GEM boundary modes out of the public biological-medium
  interface. Human and mouse workflows default explicitly to
  `normal_human_plasma` and `mouse_plasma`, respectively.
- Preserved the runnable minimal workflow at the beginning of `README.md` and
  synchronized Tutorials 1, 2, and 5, the workflow vignette, public API index,
  generated help, and medium regression tests.

# RegCompassR 2.2.0

- Made `condition_full_oof` the primary regulatory and metabolic penalty route.
  Jointly coefficient-estimable edges remain available as the common-support
  component, and the condition-unique component is reported as condition-full
  minus common support.
- Separated coefficient estimability from projection support. A candidate edge
  that is non-estimable in one or both conditions retains an unavailable
  coefficient (`NA`) but contributes an exact structural zero in each affected
  condition. Predictors equal to zero in every input cell remain represented in
  the shared candidate supergraph.
- Added the canonical Pando masks `coefficient_estimable_mask`,
  `projectable_structural_zero_mask`, and `projection_support_mask`, plus the
  public `project_condition_grn_primary_cells()` handoff.
- Removed depth matching, common-depth restriction, alpha sensitivity,
  zero-support sensitivity, and link-saturation propagation from Layer 1,
  Layer 2, persisted schemas, tests, and current documentation.
- Retained `rc_regcompass_step_target_union()` as the supported optional
  targeted-remapping pass. It scores direct KEGG/Reactome/master-Rhea-linked
  non-core reactions using the primary condition-full Layer 1 evidence and the
  exact cached Stage 5 union GEM without rebuilding the model or rerunning
  FASTCORE.
- The canonical tutorial set contains five documents: one-shot, stepwise,
  mathematical model, targeted reaction remapping, and condition comparison.
  All equations are centralized in Tutorial 3.
- Retained the existing Stage 2 graph contract: one independent graph per broad
  cell type, all conditions joint within the cell-type graph, and condition-pure
  metacells assigned after graph clustering.

# RegCompassR 1.9.3

- Added a final Pando bridge validation layer for Pando 1.2.1. RegCompass now verifies the public `infer_condition_grn()` and `condition_grn_fit()` APIs and requires the explicit `ConditionGRNFit$comparison_mask` contract.
- Condition effects can enter Layer 1 only when Pando marks the edge eligible in both the requested condition and the explicit reference condition. Complete effect tables retain all rows with `comparable_to_reference`; unsupported contrasts are excluded from active effects.
- Made `candidate_screen = "motif_domain"` the documented interaction-safe default. `condition_union` and `pooled` remain explicit marginal-screen sensitivity modes.
- Added strict ownership checks for `pando_initiate_args`, `pando_motif_args`, and `pando_infer_args`. Nested overrides of the object, assays, genome, motifs, metadata columns, GEM target genes, network name, minimum group size, error policy, and `BPPARAM` now stop before Pando is called.
- Rejected `aggregate_rna_col` and `aggregate_peaks_col` in canonical Stage 1 because Pando is fitted on paired single cells and RegCompass Stage 2 owns metacell aggregation.
- Unified Stage 1 parallel semantics. A supplied `BiocParallelParam` is forwarded to Pando; otherwise `parallel = TRUE` enables Pando native mapping, while `parallel = FALSE` is serial. `BPPARAM = TRUE` now stops explicitly. The resolved route is persisted in `step1$params$pando_parallel`.
- Pinned the merged Pando comparison-mask commit and advanced the package/result version to 1.9.3.
- Rewrote the README, one-shot and stepwise tutorials, workflow vignette, Pando mathematical contract, API index, generated help, and regression tests around the actual Pando 1.2.1 interface.
- Corrected mouse examples and help. The bundled regulatory-region objects are hg38; mouse runs must supply a species- and build-matched `GRanges` through `pando_initiate_args$regions`.

# RegCompassR 1.8.4

- Removed per-meta-module local FASTCORE from the canonical workflow. Stage 3 now produces biological meta-modules and a deduplicated merged reaction catalogue only.
- Reserved the term **union GEM** for the medium-constrained Stage 5 model. Merging meta-module reaction IDs no longer creates or names a union GEM.
- Added one global FASTCORE completion per medium-specific union GEM. The merged biological reactions are retained, and only globally required FASTCORE support reactions are added under the selected medium.
- Replaced Stage 3 `global_modules`, `global_core_reactions`, and `global_reaction_membership` outputs with `merged_modules`, `merged_core_reactions`, and `merged_reaction_membership`.
- Removed `local_completed_reaction_membership`, `local_fastcore_summary`, `local_fastcore_diagnostics`, and `local_fastcore_completion_iterations` from current workflow outputs.
- Removed the `layer1_args$local_fastcore` and `layer1_args$local_fastcore_args` interfaces. Global FASTCORE controls now live exclusively in `layer2_args$model_params`.
- Updated the canonical runner so `upstream_workers` covers GRN inference and Layer 1 only; Stage 3 no longer allocates a FASTCORE worker pool.
- Updated target-union scoring to validate anchors against the merged Stage 3 catalogue while reusing the exact cached medium-specific union GEM files.
- Synchronized README, workflow documentation, all five tutorials, the vignette, stage contracts, and generated Rd files with the global-only FASTCORE architecture.
- Added regression tests for the merged-catalogue contract, union-GEM naming, removal of local FASTCORE from Stage 3, and absence of obsolete public API names.

# RegCompassR 1.8.3

- Added a canonical two-layer worker model with `upstream_workers = 6L` for GRN/local-FASTCORE/Layer-1 tasks and `layer2_workers = 30L` for LP scoring. Setting both values to one produces a fully serial run.
- Removed `parallel_backend` from the complete-workflow interface. The package now always resolves SOCK/SnowParam on Windows and MulticoreParam on Linux/macOS automatically.
- Changed complete-run worker lifetime from a shared upstream pool to stage-scoped pools. Every package-managed pool is created for one stage, stopped on success or failure, dereferenced, and followed by `gc(full = TRUE)` before the next unrelated stage.
- Added a strict no-nested-threading contract. BLAS/OpenMP/RcppParallel and nested R worker settings are temporarily fixed at one, while Pando remains internally serial, so outer parallelism executes multiple independent single-thread analyses rather than multiplying threads inside each worker.
- Added operating-system-aware parallel configuration. Requested and actual backends, layered worker counts, one internal thread per task, OS type, stage groups, and lifecycle policy are retained in the result.
- Bundled validated Human-GEM 2.0.0 and Mouse-GEM 1.8.0 RegCompass assets under `inst/extdata/gem`. Canonical runs load them offline by default. Cache-first, explicit bundled-only, download, force-rebuild, and low-level download/update paths remain available.
- Added `rc_bundled_gem_manifest()` and exported `rc_download_species_gem()`. The installed manifest records model source, release, checksum, size, citation DOI, and CC BY 4.0 attribution.
- Added progress output and elapsed-time auditing to every public workflow stage and to the complete six-stage run. Each stage writes `step_timing.tsv`; one-shot execution writes `00_execution_timing.tsv` and stores stage and total timings in `result$timing`.
- Added an audited condition-metacell cache contract. Checkpoints are no longer reused by file existence alone: ordered cells and labels, scalable full-content RNA/ATAC fingerprints, selected PCA/LSI embedding fingerprints, the SuperCell2 label, `gamma`, seed, reductions/dimensions, and metacell thresholds must match, or the user must rebuild with `overwrite = TRUE`.
- Downstream stages now reject legacy metacell objects that lack the current condition-only label-guided construction and cache provenance instead of assigning current provenance to an unverifiable artifact.
- Kept legal minimum-version Imports for SeuratObject 4.1.4, Seurat 4.4.0, and Signac 1.11.0 while retaining the exact default versions in package Config fields; coherent v4 and v5 runtime profiles are validated separately because R dependency fields do not support profile-specific equality constraints.
- Bundled GEM loading now writes and revalidates the requested `save_rds` path, matching downloaded-model cache semantics. `rc_run_regcompass()` also preserves the prior positional location of `species` ahead of the later `progress` argument.
- Public-stage timing now records success only after the expected final RDS has been newly committed; failures during the last export/save phase are written as `status = error`.
- Hardened reaction annotation and evidence provenance: normalized GEM bounds are used, missing roles are inferred rather than forced to internal, mouse symbols retain source case, missing omnibus evidence is `unknown/unavailable`, and unavailable reaction-capacity reconstruction cannot be promoted to `RNA+ATAC` from gene-level changes alone.
- Condition plots retain full reaction annotations and evidence. Gene-associated plot collections now apply the requested condition filter to evidence selection and expose/forward `min_units` instead of using a hidden fixed value.
- Target-union scoring now validates gene and reaction selectors independently and determines target availability from the actual medium-specific cached union GEM files. Directly database-linked support reactions added during model completion remain eligible when present in all reused models, even if absent from the pre-completion membership table.
- Added a function-by-function audit of PRs #166–#171, synchronized generated help and Tutorials 3–5, and expanded regression coverage for every still-valid unresolved review finding.
- Formally documented and tested optional Harmony-based RNA geometry for Stage 2 metacells through `metacell_args$rna_reduction` and `rna_dims`; PCA remains the default, ATAC LSI remains independently selectable, and reduction names/dimensions/embedding fingerprints remain part of cache invalidation.

# RegCompassR 1.8.2

- Added `rc_regcompass_step_target_union()` for a second LP pass after the original core analysis. Selected previous core reactions are mapping anchors only. The function directly identifies non-core reactions sharing KEGG, Reactome, or master-Rhea identifiers with a selected core and scores them in the exact cached global union GEM. Same-subsystem, recursive/transitive, FASTCORE-only, generic union, and metabolite-neighbour expansion are not used; previously scored global cores are not recomputed.
- Added strict stage contracts. Layer 1 and Layer 2 now carry classes, workflow parameters, GEM fingerprints, and ordered unit identifiers; Stage 3-6 reject objects from a different GEM, workflow, or metacell order.
- Removed the retired `v170_microcompass_contract.R` compatibility override and the redundant `internal_apply.R` wrapper. Renamed the active regulatory integration helper without a historical version suffix and replaced package-version-specific algorithm labels with semantic schema identifiers.
- Updated the tutorials, vignette, README, API index, and help pages for the 1.8.2 workflow. Repetitive migration text and the obsolete architecture correction document were removed.

# RegCompassR 1.8.1

- Added formal reaction annotation to Stage 6 and condition-statistics outputs: reaction names, stoichiometry-derived formulas with metabolite names and compartments, direction-specific substrates/products, subsystems, GPR rules, participating genes, and database cross-references.
- Added condition-by-cell-type evidence provenance that distinguishes active `RNA+ATAC` support from `RNA-only`, `GPR/no-observed-RNA`, and `structural/no-GPR` reactions. `RNA+ATAC` now requires the GPR-aggregated reaction capacity calculated from integrated evidence to differ from the otherwise identical RNA-only reaction capacity; gene-level ATAC modifiers and contribution genes are reported separately.
- Added `rc_build_reaction_annotations()` and `rc_attach_reaction_annotations()` for new and previously generated results.
- Added `rc_select_gene_reactions()` and `rc_plot_condition_gene_reactions()` for selecting scored reactions by metabolic genes and generating a ranked collection of significant, biologically annotated condition boxplots.
- Added `rc_test_condition_reactions()` for same-reaction, same-direction, same-medium comparisons between conditions within each cell type. It reports Kruskal-Wallis omnibus tests, pairwise Wilcoxon tests, BH-adjusted P values, median score shifts, rank-biserial/common-language effects, and Cohen's d.
- Added `rc_plot_condition_reaction()` for multi-condition boxplots of a selected reaction target, with every metacell shown as a jittered point, Kruskal-Wallis omnibus annotation, and pairwise significance brackets based on raw or reaction-wide multiplicity-adjusted P values.
- Condition-reaction statistics explicitly distinguish within-dataset metacell significance from biological-replicate-level treatment inference and verify that target `vmax` is invariant across units before testing.
- Fixed `mouse_plasma` so it no longer inherits human HPLM concentrations or provenance. Healthy-mouse glucose (4.381 mM), lactate (3.088 mM), and glutamine (0.934 mM) define the only quantitative relative uptake caps; all other mouse components are availability-only.
- Separated the healthy-mouse quantitative reference from the broader murine plasma and tumor-interstitial-fluid availability evidence, and removed the unrelated Mouse-GEM reconstruction DOI from medium-composition provenance.
- Removed the redundant public `metacell_label_col` and stepwise `label_col` arguments. The canonical workflow now exposes its actual behavior directly: `celltype_col` is always passed to SuperCell2 before aggregation, while condition remains the only hard metacell stratum.
- Retained `label_col` only on the lower-level general-purpose `rc_make_supercell2_metacells()` builder, where it is a functional SuperCell2 option.
- Updated the README, all three tutorial levels, the workflow vignette, API index, and help pages to use the canonical interface only.
- Added a complete guide to the predefined extracellular media, including species restrictions, assumptions, and custom-medium examples.

# RegCompassR 1.7.0

- Changed the canonical metacell scope to `condition × cell type`, deliberately pooling cells from all biological samples within each condition before SuperCell2 while retaining per-metacell biological-sample composition diagnostics.
- Changed Pando inference and GRN meta-module construction to the same condition-by-cell-type scope.
- Allows Pando installed from a locally downloaded source archive when GitHub remote metadata are unavailable. Such installations continue with an explicit warning and are marked as having an unverified repository origin; explicitly conflicting remote username or repository metadata still fail.
- Uses condition-specific Pando coefficients learned from RNA+ATAC to weight accessibility-only regulatory deviations at the metacell level; metacell TF RNA is not multiplied into the modifier, reducing direct duplicate RNA weighting.
- Clarifies that coefficients estimated from the same pooled dataset are fitted parameters rather than independent validation evidence; condition-pooled outputs remain descriptive unless external fitting or cross-fitting is supplied.
- Fixed the canonical GPR calculation to a normalized, monotone Boltzmann soft-min AND, additive isozyme OR, and no promiscuity weighting. This historical rule is superseded in the development version by COMPASS `min`/`median`/`mean` aggregation.
- Replaced the previous decomposed expression-plus-confidence objective with one COMPASS-like positive cost, `1 / (1 + log2(1 + E_multiome))`.
- Restricts fixed structural penalties to exchange, demand, sink, and artificial-support reactions. Transport and cofactor reactions with GPR evidence retain the integrated multiome reaction-expression cost.
- Builds biological meta-modules only from complete-GPR core reactions, core-reaction subsystems, and reactions sharing KEGG, Reactome, or master-Rhea identifiers. Metabolite-neighbour expansion is not used; local FASTCORE is the sole mechanism for adding reactions required for flux feasibility.
- Supports both shared union meta-module GEM and shared full-GEM scoring modes with the same Layer 1 evidence, medium, target-flux fraction, and ranking outputs.
- Allows one or more biological samples per condition. Sample counts are retained as provenance and do not block the descriptive pooled-metacell workflow.
- Allows one condition. Single-condition runs return within-condition reaction priorities; multi-condition runs additionally return all pairwise descriptive priority contrasts within each cell type.
- Added explicit `reaction_ranking` output containing reaction ID, direction, medium, median minimum penalty, support score, and within-condition priority rank.
- Deleted obsolete sample-level differential/statistics code and unused pseudobulk interfaces that were incompatible with the pooled-metacell inference semantics.
- Deleted the retired strict-stratum global workflow, Q95 calibration implementation, Pando reaction-confidence implementation, Layer 2 confidence alignment functions, confidence placeholders, `penalty_weights` API, and metabolite-neighbour expansion helper and controls.
