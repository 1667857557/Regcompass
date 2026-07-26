# Stage input-output contracts

RegCompassR 1.8.8 connects stages only when classes, workflow settings, GEM provenance, metacell provenance, and scoring-unit order agree.

## Stage 1: shared-background condition sub-GRNs

Class: `regcompass_grn_step`

Required outputs in canonical mode:

```r
step1$grn_result$target_metabolic_genes
step1$grn_result$celltype_fit_status
step1$grn_result$sample_status
step1$grn_result$tf_peak_gene_candidates
step1$grn_result$tf_peak_gene_global
step1$grn_result$tf_peak_gene_condition_all
step1$grn_result$tf_peak_gene_significant
step1$grn_result$condition_target_genes
step1$grn_result$target_model_diagnostics
step1$grn_result$stability_diagnostics
step1$grn_result$normalization_policy
step1$gem_fingerprint
step1$params
```

Contract:

- one Pando candidate universe is built per cell type and identified by `edge_universe_id`;
- every condition in that cell type uses the same edge dictionary;
- `effective_estimate = global_estimate + condition_deviation`;
- condition deviations sum to zero for every edge;
- `estimate` is the stability-adjusted coefficient consumed downstream;
- an active edge satisfies the configured effect, cross-validated reliability, selection-frequency and sign-stability thresholds;
- `padj` is `NA` in multitask mode and `evidence_type` records the stability-selection policy;
- every `condition × cell type` group has one `sample_status` row, even when it has no active edges.

The default structural regions remain species specific:

```text
human = union(Pando::phastConsElements20Mammals.UCSC.hg38,
              Pando::SCREEN.ccRE.UCSC.hg38)
mouse = Pando::phastConsElements20Mammals.UCSC.hg38
```

## Stage 2: metacells

Class: `regcompass_metacell_step`

Required outputs:

```r
step2$pooled$metacell_meta
step2$pooled$membership
step2$metacell_object
step2$params
```

The merged metacell object and metadata must contain the same ordered units. Condition remains the hard pooling stratum; cell type is used as the SuperCell2 label and audited after aggregation. Biological sample composition remains provenance.

## Stage 3: condition biological meta-modules

Class: `regcompass_meta_module_step`

Required outputs:

```r
step3$condition_modules$supported_metabolic_genes
step3$condition_modules$core_gene_reaction
step3$condition_modules$reaction_membership
step3$condition_modules$meta_module_summary
step3$merged_modules$merged_core_reactions
step3$merged_modules$merged_reaction_membership
step3$merged_modules$source_edge_universe_ids
step3$group_coverage
```

Contract:

- `supported_metabolic_genes` contains condition sub-GRN targets with active stable edges;
- positive and negative stable edges both establish regulated-gene membership;
- a reaction is core only when one complete GPR branch is contained in the condition target set;
- partial GPR complexes cannot anchor expansion;
- expansion is one ordered pass: core subsystem, direct KEGG/Reactome equivalence, then direct master-Rhea equivalence;
- `merged_core_reactions` and `merged_reaction_membership` are deduplicated reaction catalogues;
- the merged object is not a GEM, has `is_gem = FALSE`, and has `fastcore_applied = FALSE`;
- Stage 3 does not apply medium bounds or run FASTCORE.

The only optional Stage 3 parameter is:

```r
meta_module_args = list(
  subsystem_table = custom_subsystem_table
)
```

## Stage 4: Layer 1

Class: `regcompass_layer1_step`

Required outputs:

```r
step4$gene_support_rna
step4$gene_regulatory_modifier
step4$gene_support_multiome
step4$reaction_expression
step4$metacell_meta
step4$capacity_params$and_method
step4$workflow_params
step4$gem_fingerprint
```

The regulatory modifier:

- uses Stage 1 stable condition coefficients;
- projects only metacell ATAC deviations;
- signed-sums TFs sharing one measured peak and target;
- uses one target denominator shared across conditions;
- uses joint target `cv_rsq` as reliability;
- gives a zero modifier to genes without active condition edges.

The zero modifier exactly recovers RNA-only support under the log-odds update. `capacity_params$and_method` must be one of `min`, `median`, or `mean`; the default is `min`.

## Stage 5: Layer 2

Class: `regcompass_layer2_step`

For `model_mode = "meta_module_gem"`, Stage 5 requires:

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

For each medium, Stage 5 constructs one final union GEM and performs one global FASTCORE completion. All conditions and metacells in that medium must resolve to the same model file and therefore the same reaction IDs, stoichiometric matrix, lower bounds and upper bounds.

Different media may produce different support-completed structures. Within a medium, structural variation across conditions is prohibited.

FASTCORE controls are supplied only through:

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

`completion_time_limit` limits union-GEM construction, not scoring LPs.

## Stage 6: results

The final result contains:

```r
result$version
result$grn_mode
result$grn
result$condition_grn_meta_modules
result$merged_grn_meta_modules
result$microcompass
result$reaction_ranking
result$condition_summary
result$condition_contrast
```

`result$grn` preserves the full Stage 1 candidate/global/condition/stability contract. `merged_grn_meta_modules` remains the Stage 3 catalogue. `microcompass$model_cache_summary` identifies the final Stage 5 union GEMs.

## Target-union restart contract

`rc_regcompass_step_target_union()` requires the original Stage 3 merged catalogue, original Layer 1 matrix, completed `meta_module_gem` Stage 5 object, and accessible model files with matching checksums and medium identifiers. It reuses the exact cached final union GEM and does not rerun FASTCORE or change medium bounds.

See [multitask GRN mathematics and object contracts](multitask-shared-grn.md).
