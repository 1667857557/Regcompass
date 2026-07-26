# RegCompassR 1.8.8

- Added the canonical `multitask_shared_backbone` GRN mode. Pando supplies one validated, condition-agnostic structural TF–peak–target candidate universe per cell type, and RegCompass jointly estimates a global backbone plus symmetric zero-sum condition deviations.
- Added condition-balanced elastic-net fitting. Every condition contributes the same total loss weight, and `alpha < 1` is required so the centred condition-deviation parameterisation has a unique ridge-regularised solution.
- Replaced repeated subsampling with full-size condition-stratified nonparametric bootstrap. Each condition is resampled with replacement at its original cell count, then the target and TF-by-ATAC predictors are re-centred within the bootstrap condition before fitting at the full-data selected lambda.
- Defined edge stability by bootstrap selection frequency and conditional sign stability. The Layer 1 coefficient is the full-data condition effect multiplied by both reliability terms; active-edge membership remains threshold based.
- Removed the biological-sample column from the canonical public workflow. Stage 1 residualisation, cross-validation, bootstrap, Stage 2 metacells, tutorials, Rd files, and result provenance are condition based only.
- Added explicit Stage 1 output contracts for structural candidates, global coefficients, all condition coefficients, bootstrap-active condition edges, condition target genes, model diagnostics, and bootstrap diagnostics.
- Replaced condition-wise coefficient normalisation in the Layer 1 regulatory modifier. TFs sharing one measured peak are signed-summed, and one target-specific denominator is shared across conditions, preserving differences in total regulatory strength.
- Preserved exact RNA-only fallback: a gene without a bootstrap-active edge in one condition has regulatory modifier zero and therefore retains its unmodified RNA support.
- Kept complete-GPR condition core construction and the one-pass subsystem/KEGG–Reactome/master-Rhea expansion. Internal and exported module identifiers now use `group_id` for the actual `condition × cell type` analysis unit.
- Kept one medium-specific union GEM per medium. All conditions and metacells reuse identical reaction IDs, stoichiometry, lower bounds, and upper bounds; only evidence-derived penalties differ.
- Added mathematical and regression tests for symmetric coding, zero-sum deviations, equal condition weights, full-size bootstrap sampling with replacement, bootstrap re-centring, no-sample public signatures, RNA-only fallback, condition-specific complete-GPR cores, and merged-reaction provenance.
- Added `glmnet` as a direct dependency and raised the required Pando fork to version 1.1.3 for the validated version-2 structural design contract.

# RegCompassR development

- Added `rc_report_condition_directions()` as a final reporting layer that retains forward/reverse LP targets, diagnoses numerically indistinguishable directions, and derives non-additive `any_direction_support` and `directional_balance` summaries. The latter is explicitly support asymmetry rather than net flux.
- Added explicit Seurat runtime profiles: the pinned SeuratObject 4.1.4 / Seurat 4.4.0 / Signac 1.11.0 stack remains the canonical default, while coherent SeuratObject/Seurat 5.x stacks with Signac >=1.12.0,<2 are accepted as a compatibility profile.
- Added assay-class-aware matrix access for v3 `Assay`, Signac 1.x `ChromatinAssay`, and joinable Seurat v5 `Assay5` objects. Split `counts.*` and `data.*` layers are joined on a working copy and recorded as provenance; ambiguous layer layouts and Signac 2.x `ChromatinAssay5` stop explicitly.
- Redefined Stage 3 around the actual analysis target: for each `condition × cell type`, GEM metabolic genes with at least one significant Pando TF–peak–gene coefficient form one supported gene set, and reactions become core only when a complete GPR branch is contained in that set.
- Removed the retired shared-TF target projection, signed target-target component construction, top-k pruning, TF-Jaccard filtering, and per-TF target truncation code. The associated functions, output fields, parameters, tests, and documentation were deleted rather than retained as compatibility wrappers.
- Added `supported_metabolic_genes`, a condition-by-cell-type gene evidence table reporting significant edge, TF, region, adjusted-P-value, coefficient, model-R², and positive/negative edge summaries.
- Made Pando's bundled `motifs` data object the canonical default when `pfm` is omitted. Explicit user-supplied motif collections remain supported.
- Added species-specific Pando region defaults. Unless overridden by `pando_initiate_args$regions`, human Stage 1 loads `phastConsElements20Mammals.UCSC.hg38` and `SCREEN.ccRE.UCSC.hg38` from Pando and uses their union, while mouse Stage 1 uses only `phastConsElements20Mammals.UCSC.hg38`.
- Replaced the former Boltzmann GPR-AND calculation with the three COMPASS aggregation functions `min`, `median`, and `mean`. The canonical default is `min`; the Boltzmann helper, `tau` parameter, associated tests, and compatibility paths were deleted.
- Fixed Stage 3 expansion to one ordered pass: core-subsystem reactions, direct KEGG/Reactome reaction equivalents, then direct master-Rhea reaction equivalents. Removed `expansion_mode`, `max_iterations`, fixed-point recursion, iteration output fields, and any one-hop/stoichiometric-neighbour expansion interface.
- Split canonical configuration into `meta_module_args` for an optional custom subsystem table and `layer1_args` for Stage 4 integrated-evidence parameters.
- Reordered public runner arguments by processing sequence: shared inputs, Stage 1 Pando, Stage 2 metacells, Stage 3 meta-modules, Stage 4 Layer 1, Stage 5 Layer 2, and execution controls.
- Removed scoring `time_limit` from directional LP and second-pass APIs. `layer2_args$model_params$completion_time_limit` remains exclusively for FASTCORE construction of the medium-specific union GEM.
- Expanded the main and stepwise tutorials with Pando evidence filters and complete metacell geometry/reproducibility settings.
- Updated the README, stepwise and one-shot tutorials, workflow vignette, mathematical workflow, stage contracts, generated help, and regression tests.

# RegCompassR 1.8.4

- Removed per-meta-module local FASTCORE from the canonical workflow. Stage 3 now produces biological meta-modules and a deduplicated merged reaction catalogue only.
- Reserved the term **union GEM** for the medium-constrained Stage 5 model. Merging meta-module reaction IDs no longer creates or names a union GEM.
- Added one global FASTCORE completion per medium-specific union GEM. The merged biological reactions are retained, and only globally required FASTCORE support reactions are added under the selected medium.
- Replaced Stage 3 `global_modules`, `global_core_reactions`, and `global_reaction_membership` outputs with `merged_modules`, `merged_core_reactions`, and `merged_reaction_membership`.
- Removed local FASTCORE output fields and controls from Stage 3.
- Updated target-union scoring to validate anchors against the merged Stage 3 catalogue while reusing the exact cached medium-specific union GEM files.

# RegCompassR 1.8.3

- Added a canonical two-layer worker model with stage-scoped worker pools and strict no-nested-threading behavior.
- Bundled validated Human-GEM 2.0.0 and Mouse-GEM 1.8.0 assets for offline use.
- Added progress, timing, cache-contract, Seurat compatibility, and reaction-provenance auditing.

# RegCompassR 1.8.2

- Added targeted second-pass reaction scoring in the exact cached final union GEM and strict stage contracts.

# RegCompassR 1.8.1

- Added reaction annotation, condition statistics, evidence provenance, predefined media documentation, and condition-level metacell workflow support.

# RegCompassR 1.7.0

- Introduced condition-pooled multiome GRN and metacell analysis, integrated RNA+ATAC reaction support, complete-GPR meta-modules, shared structural scoring, and directional reaction ranking.
