# RegCompassR 1.8.8

- Added the default `multitask_shared_backbone` GRN mode. Pando now constructs one structural TF–peak–metabolic-target candidate universe per cell type across all conditions, and RegCompass fits condition effects on identical design columns.
- Added a condition-balanced elastic-net model with a common global block and more strongly penalized condition-specific blocks. The exported reference-free decomposition is `effective = global + condition_deviation`, with global defined as the condition mean and deviations constrained to sum to zero after canonicalization.
- Added condition-stratified cross-validation and stability subsampling. Condition edges report selection frequency, sign stability, stability weight, active-edge state, effective direction and stable sign-reversal flags instead of treating regularized coefficients as classical adjusted-p-value tests.
- Added explicit `tf_peak_gene_candidates`, `tf_peak_gene_global`, `tf_peak_gene_condition_all`, `condition_target_genes`, `celltype_fit_status`, and shared Pando `design_id` outputs.
- Preserved the complete TF–peak–gene chain through explicit source RNA/ATAC feature identifiers exported by Pando 1.1.2.
- Defined condition-regulated metabolic genes from active positive or negative edges. Stage 3 retains the strict complete-GPR rule: a reaction becomes a condition core only when at least one complete GPR AND branch is contained in that condition/cell-type target set.
- Replaced condition-local multitask edge normalization in Layer 1 with a target/cell-type scale shared across conditions. Effective coefficients are multiplied by stability, a condition-balanced TF reference and the original interaction scale; TF edges sharing one ATAC peak are signed-summed before accessibility projection.
- Preserved `legacy_condition_pando` as an explicit compatibility mode, including the former adjusted-p-value filtering and condition-local Layer 1 weight normalization.
- Retained the Stage 5 structural contract: each medium produces one shared union GEM with one global FASTCORE completion, and the exact stoichiometric matrix, bounds and target catalogue are reused for every condition and metacell.
- Added mathematical-invariant, schema, compatibility and workflow regression tests; added `glmnet` as a direct dependency; updated the README and shared-backbone mathematical documentation.

# RegCompassR development before 1.8.8

- Added `rc_report_condition_directions()` as a final reporting layer that retains forward/reverse LP targets, diagnoses numerically indistinguishable directions, and derives non-additive `any_direction_support` and `directional_balance` summaries. The latter is explicitly support asymmetry rather than net flux.
- Added explicit Seurat runtime profiles: the pinned SeuratObject 4.1.4 / Seurat 4.4.0 / Signac 1.11.0 stack remains the canonical default, while coherent SeuratObject/Seurat 5.x stacks with Signac >=1.12.0,<2 are accepted as a compatibility profile.
- Added assay-class-aware matrix access for v3 `Assay`, Signac 1.x `ChromatinAssay`, and joinable Seurat v5 `Assay5` objects. Split `counts.*` and `data.*` layers are joined on a working copy and recorded as provenance; ambiguous layer layouts and Signac 2.x `ChromatinAssay5` stop explicitly.
- Redefined Stage 3 around the actual analysis target: for each `condition × cell type`, GEM metabolic genes with regulatory TF–peak–gene evidence form one supported gene set, and reactions become core only when a complete GPR branch is contained in that set.
- Removed the retired shared-TF target projection, signed target-target component construction, top-k pruning, TF-Jaccard filtering, and per-TF target truncation code.
- Added species-specific Pando motif and regulatory-region defaults.
- Replaced the former Boltzmann GPR-AND calculation with COMPASS-compatible `min`, `median`, and `mean`; the canonical default is `min`.
- Fixed Stage 3 expansion to one ordered pass: core-subsystem reactions, direct KEGG/Reactome equivalents, then direct master-Rhea equivalents.
- Removed scoring `time_limit` from directional LP and second-pass APIs. `layer2_args$model_params$completion_time_limit` remains exclusively for FASTCORE construction.

# RegCompassR 1.8.4

- Removed per-meta-module local FASTCORE from the canonical workflow. Stage 3 now produces biological meta-modules and a deduplicated merged reaction catalogue only.
- Reserved the term **union GEM** for the medium-constrained Stage 5 model. Merging meta-module reaction IDs no longer creates or names a union GEM.
- Added one global FASTCORE completion per medium-specific union GEM. The merged biological reactions are retained, and only globally required FASTCORE support reactions are added under the selected medium.
- Replaced Stage 3 `global_modules`, `global_core_reactions`, and `global_reaction_membership` outputs with `merged_modules`, `merged_core_reactions`, and `merged_reaction_membership`.
- Removed local FASTCORE outputs and interfaces from Stage 3 and Layer 1.
- Updated target-union scoring to validate anchors against the merged Stage 3 catalogue while reusing the exact cached medium-specific union GEM files.

# RegCompassR 1.8.3

- Added a canonical two-layer worker model with `upstream_workers = 6L` for GRN and Layer-1 tasks and `layer2_workers = 30L` for LP scoring.
- Added stage-scoped worker pools, a strict no-nested-threading contract, operating-system-aware parallel configuration, progress output, elapsed-time auditing, and audited metacell cache fingerprints.
- Bundled validated Human-GEM 2.0.0 and Mouse-GEM 1.8.0 assets for offline use.
- Hardened reaction annotation, evidence provenance, medium handling and target-union scoring.

# RegCompassR 1.8.2

- Added `rc_regcompass_step_target_union()` for a second LP pass using direct KEGG, Reactome, or master-Rhea links while reusing the exact cached union GEM.
- Added strict stage classes, workflow parameters, GEM fingerprints and unit-order validation.
- Updated tutorials, vignette, README, API index and help pages.

# RegCompassR 1.8.1

- Added reaction annotations, evidence provenance, gene–reaction selection and plotting, condition reaction statistics and condition reaction plots.
- Added medium presets and removed redundant public metacell label arguments.

# RegCompassR 1.7.0

- Introduced condition-pooled metacells, condition-specific regulatory evidence, COMPASS-like reaction penalties and shared union/full-GEM scoring modes.
- Clarified that coefficients learned from the same pooled data are fitted parameters rather than independent validation evidence.
