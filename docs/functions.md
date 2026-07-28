# Public functions and API contract in RegCompassR 1.9.1

RegCompassR exposes 23 functions. The canonical analysis is available as a
one-shot call or as six validated, restartable stages. Pando is the only GRN
estimator; RegCompass consumes its versioned `ConditionGRNFit` without
refitting condition networks.

## Model, medium, and execution setup

| Function | Purpose | Main return value |
|---|---|---|
| `rc_prepare_gem()` | Load a pinned Human-GEM or Mouse-GEM from cache/bundled assets, or explicitly download and rebuild it. | Validated RegCompass GEM |
| `rc_prepare_human2_gem()` | Human-GEM 2.0.0 convenience entry point. | Validated human GEM |
| `rc_prepare_mouse_gem()` | Mouse-GEM 1.8.0 convenience entry point; mouse symbols are retained. | Validated mouse GEM |
| `rc_bundled_gem_manifest()` | Inspect bundled release, checksum, size, source, citation, and license metadata. | Manifest data frame |
| `rc_download_species_gem()` | Download and parse an official upstream model for asset rebuilding or a newer pinned release. | Parsed model tables and download diagnostics |
| `rc_make_medium_scenarios()` | Build a condition-invariant extracellular medium table. | Reaction-level medium constraints |
| `rc_parallel_config()` | Report platform-aware backend and worker resolution without starting workers. | Execution configuration list |

## Complete workflows

| Function | Purpose |
|---|---|
| `rc_run_regcompass_one_shot()` | Prepare the species GEM and medium when omitted, then run the canonical workflow. |
| `rc_run_regcompass()` | Run Stages 1–6 with explicit GEM, medium, worker budgets, and argument bundles. |

`rc_run_regcompass()` accepts:

- `pando_args` for Stage 1;
- `metacell_args` for Stage 2;
- `meta_module_args` for Stage 3;
- `layer1_args` for Stage 4;
- `layer2_args` and `model_mode` for Stage 5;
- `upstream_workers` and `layer2_workers` for execution.

`sample_col` remains in the complete-run and Stage 2 signatures only for
backward call compatibility. It is not used for canonical grouping, sample
selection, balancing, weighting, stability selection, or refitting.

## Restartable stages

| Stage | Function | Input → output |
|---:|---|---|
| 1 | `rc_regcompass_step_grn()` | Paired RNA+ATAC cells and GEM → one `ConditionGRNFit` per cell type plus condition coefficient/effect tables |
| 2 | `rc_regcompass_step_metacells()` | Paired RNA+ATAC cells → condition-pooled, cell-type-guided SuperCell2 metacells |
| 3 | `rc_regcompass_step_meta_modules()` | Active condition target genes → complete-GPR cores and a deduplicated biological reaction catalogue |
| 4 | `rc_regcompass_step_layer1()` | Metacells, Pando transforms/effects, and GPRs → RNA-only and RNA+ATAC reaction-support matrices |
| 5 | `rc_regcompass_step_layer2()` | Layer 1 support and shared medium → cached shared structural GEMs and directional LP scores |
| 6 | `rc_regcompass_step_results()` | All validated stages → annotated `regcompass_condition_grn_fit_v2` result |
| optional | `rc_regcompass_step_target_union()` | Any GEM reaction-ID anchor or an original-core gene selector → directly linked non-core targets scored in the exact cached Stage 5 model |

Each stage validates its input class, workflow settings, GEM fingerprint, and
ordered units. Each writes `step_timing.tsv` and a classed RDS checkpoint.

## Stage 1: authoritative condition-GRN contract

The required defaults are:

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
    reference_condition = "Control",
    scale = TRUE
  )
)
```

Conditions within one cell type share the complete edge dictionary,
eligibility mask, pooled final `TF RNA × peak ATAC` transform, target scale,
lambda path, and selected lambda. Condition coefficients are estimated
independently. The exported effect is
`beta_condition - beta_reference`; the Universal coefficient mean is retained
only for Pando visualization compatibility.

Current Stage 1 outputs include:

```r
step1$grn_result$condition_grn_fits
step1$grn_result$condition_fit_status
step1$grn_result$tf_peak_gene_condition
step1$grn_result$tf_peak_gene_condition_effect
step1$grn_result$normalization_policy
```

`sample_status`, `tf_peak_gene_all`, and `tf_peak_gene_significant` remain
lossless compatibility aliases. New code should use the current fields above.
Stage 1 active-edge selection uses coefficient magnitude and finite
target-model R²; it does not use coefficient-level adjusted P values.

## Stage 2: geometry and cache contract

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

Cells are stratified by condition only. Cell type is supplied as a
construction label and audited by dominant membership after aggregation. The
selected reductions, dimensions, embedding fingerprints, ordered cells,
assays, seed, gamma, and thresholds are part of the cache contract.

## Stages 3–5

Stage 3 maps active condition target genes to complete-GPR cores and performs
one ordered annotation pass:

```text
core subsystem
→ direct KEGG/Reactome reaction equivalence
→ direct master-Rhea equivalence
```

Stage 3 does not run FASTCORE and does not construct a GEM.

Stage 4 reconstructs the standardized TF×ATAC predictor from the Pando
transform and projects the explicit reference contrast. It does not refit a
GRN and does not divide effects by their absolute sum. `gpr_and_method` accepts
`"min"`, `"median"`, or `"mean"` and defaults to `"min"`; isozyme OR branches
remain additive.

With `model_mode = "meta_module_gem"`, Stage 5 applies each medium, runs the
single global FASTCORE completion, caches one union GEM, and reuses that exact
model for every condition and metacell in the medium. With
`model_mode = "full_gem"`, the validated complete GEM is reused without
union-model reconstruction.

## Annotation, statistics, and plots

| Function | Purpose |
|---|---|
| `rc_build_reaction_annotations()` | Build names, directional formulas, GPR genes, database identifiers, and RNA/RNA+ATAC evidence provenance. |
| `rc_attach_reaction_annotations()` | Add the current annotation contract to an existing RegCompass result. |
| `rc_select_gene_reactions()` | Select scored reactions containing any or all requested metabolic genes. |
| `rc_test_condition_reactions()` | Compare the same reaction, direction, medium, and cell type across conditions. |
| `rc_report_condition_directions()` | Preserve forward/reverse results and derive non-additive best-direction and direction-balance summaries. |
| `rc_plot_condition_reaction()` | Plot one fixed reaction-direction target across conditions. |
| `rc_plot_condition_gene_reactions()` | Select reaction targets by metabolic genes and return annotated condition plots. |

Condition tests use metacells as descriptive pseudo-observations. They do not
create biological replicates or estimate sample-level treatment effects.
Forward and reverse objectives are separate LP targets and are never added or
interpreted as measured net flux.

## Documentation

- [Tutorial index](run-modes-and-stepwise-workflow.md)
- [Level 1: one-shot workflow](tutorial-01-quick-start.md)
- [Level 2: stepwise workflow](tutorial-02-stepwise-audit.md)
- [Level 3: restart and diagnostics](tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition comparison](tutorial-05-condition-differential-analysis.md)
- [Condition-comparable GRN and penalty contract](condition-comparable-grn.md)
- [Stage input-output contracts](stage-interface-contracts.md)
