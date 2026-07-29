# Stage input-output contracts

RegCompassR 1.9.1 connects stages only when classes, workflow settings, GEM provenance, metacell construction provenance, and scoring-unit order agree.

## Stage 1: Pando evidence

Class: `regcompass_grn_step`

Required outputs:

```r
step1$grn_result$target_metabolic_genes
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$normalization_policy$pando_motifs
step1$grn_result$normalization_policy$pando_regions
step1$grn_result$normalization_policy$pando_evidence_filters
step1$gem_fingerprint
step1$params
```

Every condition represented within a cell-type `ConditionGRNFit` must complete
successfully. A successful condition layer may legitimately have zero active
coefficients; only groups with active supported target genes can contribute
complete-GPR cores. `target_metabolic_genes` is the intersection of GEM GPR
genes and RNA-assay row names.

The complete fit contract must use schema
`pando_condition_grn_fit_v2`, the
`shared_baseline_condition_sparse_elastic_net` engine, pooled final-edge and target
standardization, one explicit `reference_condition`, aligned `beta`,
`contrast`, and `eligibility_mask` matrices, and stored predictor/response
transforms. The contrast must equal
`beta_condition - beta_reference`.

`sample_status`, `tf_peak_gene_all`, and `tf_peak_gene_significant` are retained
as compatibility aliases. New code should use the current fields listed above.

When `pfm` is omitted, `pando_motifs` records `Pando::motifs`, loaded with `data("motifs", package = "Pando")`. Without an explicit `pando_initiate_args$regions`, the region contract is species-specific:

```text
human = union(Pando::phastConsElements20Mammals.UCSC.hg38,
              Pando::SCREEN.ccRE.UCSC.hg38)
mouse = Pando::phastConsElements20Mammals.UCSC.hg38
```

`step1$params$species` records the resolved species, and `pando_regions` records the applied default or `user_supplied` policy.

## Stage 2: metacells

Class: `regcompass_metacell_step`

Required outputs:

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$metacell_object
step2$params
```

The merged metacell object and metadata must contain the same ordered units. Reduction names, dimensions, cell labels, assay fingerprints, and embedding fingerprints are part of the cache contract.

Cells are grouped by condition only. Cell type is a construction label followed
by dominant-membership auditing. User sample metadata do not enter selection,
weighting, grouping, stability selection, or model refitting.

## Stage 3: biological meta-modules

Class: `regcompass_meta_module_step`

Required outputs:

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$condition_modules$meta_module_summary
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$group_coverage
```

Contract:

- `supported_metabolic_genes` contains one row per condition, cell type, and GEM target gene with at least one active Pando TF–peak–gene row;
- positive and negative Pando coefficients both count as regulatory evidence;
- all supported genes in one `condition × cell type` form one GPR-evaluation set;
- Stage 3 performs no shared-TF target projection, top-k graph pruning, or connected-component analysis;
- `core_gene_reaction` marks a reaction as core only when one complete GPR branch is contained in the supported gene set;
- expansion is exactly one ordered pass: core subsystem, direct KEGG/Reactome equivalence, then direct master-Rhea equivalence;
- a reaction added at the master-Rhea step does not trigger another subsystem or KEGG/Reactome pass;
- `merged_core_reactions` contains deduplicated complete-GPR core reactions;
- `merged_reaction_membership` contains deduplicated biological reactions only;
- `merged_modules$is_gem` is `FALSE`;
- `merged_modules$fastcore_applied` is `FALSE`;
- Stage 3 does not apply medium constraints or run FASTCORE;
- the merged object is a catalogue and must not be described as a union GEM.

The only optional Stage 3 parameter is a custom subsystem table:

```r
meta_module_args = list(
  subsystem_table = custom_subsystem_table
)
```

Omitting `meta_module_args` uses the GEM's subsystem annotations.

## Stage 4: Layer 1

Class: `regcompass_layer1_step`

Required outputs:

```r
step4$reaction_expression
step4$metacell_meta
step4$capacity_params$and_method
step4$workflow_params
step4$gem_fingerprint
```

The reaction-expression matrix must contain every merged core reaction and the same ordered metacells represented by Stage 2.

The gene-level modifier must be reconstructed from the stored Pando transform
of metacell `TF RNA × peak ATAC` and the explicit reference contrast. Layer 1
must not refit the GRN, use the Universal row mean as a baseline, or normalize
condition effects by their absolute sum.

`capacity_params$and_method` must be one of:

```r
c("min", "median", "mean")
```

The canonical default is `"min"`.

## Stage 5: Layer 2

Class: `regcompass_layer2_step`

For `model_mode = "meta_module_gem"`, required inputs are:

```r
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
```

Required outputs include:

```r
step5$model_cache_summary
step5$source_core_reactions
step5$source_merged_reaction_membership
step5$score
step5$penalty
step5$vmax
```

Each `model_cache_summary` row identifies one final medium-specific union GEM and records at least:

```r
medium_scenario
file
file_checksum
build_strategy
completion_stage
```

The cached model must record:

```r
model$is_union_gem
model$union_gem_medium_scenario
model$build_params$completion_stage
model$reaction_meta$global_fastcore_support
```

The required build metadata are:

```text
build_strategy = medium_specific_union_gem
completion_stage = single_global_fastcore_after_meta_module_merge
```

All conditions and metacells within the same medium must resolve to the same model file. Different media may resolve to different union-GEM structures because their global FASTCORE support sets may differ.

## Stage 6: results

The final result contains:

```r
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$grn$condition_grn_fits
result$grn$condition_fit_status
result$grn$tf_peak_gene_condition
result$grn$tf_peak_gene_condition_effect
result$merged_grn_meta_modules
result$microcompass
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

`merged_grn_meta_modules` is the Stage 3 catalogue. `microcompass$model_cache_summary` identifies the final Stage 5 union GEMs.

## Global FASTCORE configuration

The only structural completion controls are supplied at Stage 5:

```r
layer2_args = list(
  model_params = list(
    completion_time_limit = 600,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE
  )
)
```

## Target-union restart contract

`rc_regcompass_step_target_union()` requires:

- the original Stage 3 merged catalogue for anchor provenance;
- the original Layer 1 matrix;
- the completed `meta_module_gem` Stage 5 object;
- accessible final Stage 5 union-GEM files;
- matching model-file checksums and medium identifiers.

Selected genes resolve only original complete-GPR cores. Reaction-ID anchors
may be any valid reaction in the supplied GEM. Both determine direct mapping
anchors only. Target availability and all LP calculations are evaluated in the
exact cached final union GEMs. The second pass does not rebuild a GEM, change
medium bounds, or rerun FASTCORE.
