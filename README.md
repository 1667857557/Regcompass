# RegCompassR

RegCompassR 2.0.0 implements a condition-comparable regulatory–metabolic
workflow for paired single-cell RNA+ATAC data.

## Canonical architecture

```text
cells of one cell type across conditions
→ one Pando motif/domain TF–peak–GEM-gene dictionary
→ pooled TF-RNA × peak-ATAC predictor transforms
→ condition-sparse selection and common-metric refit
→ sample-blocked OOF reliability
→ active metabolic targets and complete-GPR core reactions
→ one ordered subsystem / KEGG–Reactome / master-Rhea expansion
→ SuperCell metacells within condition × broad-cell-type strata
→ single-cell TF×ATAC projection followed by membership aggregation
→ RNA+ATAC reaction support and penalties
→ one shared medium-specific GEM
→ directional COMPASS-like LP scoring
```

Pando is the sole GRN estimator. RegCompass consumes Pando 1.4.0's versioned
`ConditionGRNFit v4` without refitting its coefficient matrices. For the main
penalty path, the target-gene projection for cell `i`, gene `g`, and condition
`c` is

```text
G[i, g, c] = sum_e z[i, e] * beta_condition[e, c]
```

Two-condition quantitative comparisons use the pairwise intersection of
estimable edges. Reference contrasts remain available for interpretation but
do not enter the main metabolic penalty. Pando calculates TF×ATAC and applies
its pooled transform in single cells; RegCompass then averages the completed
gene projections using SuperCell's exact
`misc$membership_table(cell_id, metacell_id)`.

## Installation

### Default validated profile: Seurat v4

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github("1667857557/SuperCell_Seurat_V4@Supercell2")
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

### Optional compatible profile: Seurat v5

A coherent SeuratObject/Seurat 5.x stack with Signac 1.12–1.x is accepted. For
v3-style assay behavior under Seurat 5:

```r
options(Seurat.object.assay.version = "v3")
```

Signac 2.x `ChromatinAssay5` is not yet supported. See
[Seurat compatibility](docs/seurat-compatibility.md) for the validated version
matrix and assay constraints.

## Minimal human one-shot run

```r
library(RegCompassR)
library(BSgenome.Hsapiens.UCSC.hg38)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem,
  scenario = "high_glucose",
  species = "human"
)

result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "Group",
  celltype_col = "cell_type",
  sample_col = "sample_id",

  pando_args = list(
    min_cells = 300L,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    pando_infer_args = list(
      method = "shared_baseline_condition_sparse",
      candidate_screen = "motif_domain",
      tf_cor = 0.1,
      peak_cor = 0,
      alpha = 0.5,
      condition_mix = 0.5,
      condition_weight = "equal",
      cv_block_col = "sample_id",
      reference_condition = "Control",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE
    )
  ),

  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 300L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  layer1_args = list(
    projection_component = "condition",
    comparison_support = "auto",
    projection_mode = "metacell_specific",
    regulatory_reliability = "sample_blocked_oof",
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),
  medium_scenarios = medium_scenarios,
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  ),
  upstream_workers = 6L,
  layer2_workers = 30L
)
```

Do not put `parallel` inside `pando_infer_args` in the canonical examples.
The one-shot workflow controls Stage 1 through `upstream_workers`, keeps Pando's
native nested loop serial, and supplies a stage-scoped `BiocParallelParam`.

## Pando argument routing

RegCompass exposes three nested bundles and validates their ownership:

| RegCompass field | Forwarded to | Typical user controls |
|---|---|---|
| `pando_initiate_args` | `Pando::initiate_grn()` | `regions`, `exclude_exons` |
| `pando_motif_args` | `Pando::find_motifs()` | motif matching controls |
| `pando_infer_args` | `Pando::infer_condition_grn()` | `candidate_screen`, `reference_condition`, elastic-net/CV controls |

RegCompass owns and rejects nested overrides of the Seurat object, assay names,
genome, motif object, condition/cell-type columns, GEM target genes, network
name, minimum condition size, group-error policy, and `BPPARAM`.
`aggregate_rna_col` and `aggregate_peaks_col` are also rejected: canonical Stage
1 fits paired single cells, while Stage 2 owns metacell aggregation.

The interaction-safe default is `candidate_screen = "motif_domain"`. It retains
motif/domain-supported candidates and lets elastic-net regularization select
edges based on the fitted `TF RNA × peak ATAC` predictor. The optional
`pooled_within_condition` applies marginal TF-target and peak-target screening
and should be treated as a sensitivity analysis. Because that screen uses the
response before cross-validation, its OOF score is not used as confirmatory
Layer 1 reliability; RegCompass sets `q = 0` for that sensitivity path.

### Conditions with one biological sample

`sample_col` is required for provenance and Pando cross-validation, but a
condition may contain only one sample. Pando then uses cell-level folds only
to select lambda and records sample-blocked OOF performance as unavailable.
RegCompass retains the fitted GRN as exploratory edge evidence, sets its
regulatory reliability to zero, and falls back to RNA-only support in Layer 1.
It never treats cells, mixed-sample metacells, or multiple metacells from one
sample as biological replicates.

Conditions with at least two samples continue to use strict sample-blocked OOF
validation. SuperCell construction is unchanged in both cases: the only hard
strata are condition × broad cell type (`condition_col × celltype_col`);
sample IDs are retained solely for composition diagnostics and GRN validation.

## Human and mouse regulatory regions

When `pfm` is omitted, RegCompass loads Pando's bundled `motifs` object. Human
analyses also default to the union of Pando's hg38 phastCons and SCREEN ccRE
regions.

Pando does not bundle a mouse-coordinate replacement for that hg38 region set.
Mouse analyses must pass a build-matched `GRanges` object whose coordinates
match both the ATAC assay and `genome`:

```r
library(BSgenome.Mmusculus.UCSC.mm10)

mouse_regions <- readRDS("mm10_regulatory_regions.rds")
stopifnot(methods::is(mouse_regions, "GenomicRanges"))

mouse_result <- rc_run_regcompass_one_shot(
  object = A_mouse,
  outdir = "RegCompass_mouse",
  genome = BSgenome.Mmusculus.UCSC.mm10,
  species = "mouse",
  pando_args = list(
    pando_initiate_args = list(regions = mouse_regions),
    pando_infer_args = list(
      candidate_screen = "motif_domain",
      reference_condition = "Control"
    )
  )
)
```

RegCompass stops rather than silently applying hg38 regions to mouse peaks.

## Inspect Stage 1

```r
result$grn$condition_grn_fits
result$grn$condition_fit_status
result$grn$tf_peak_gene_condition
result$grn$tf_peak_gene_condition_effect
result$grn$tf_peak_gene_condition_effect_all[, c(
  "edge_id", "condition", "cell_type", "condition_effect",
  "comparable_to_reference"
)]
result$grn$normalization_policy[c(
  "pando_candidate_screen", "comparison_support",
  "reference_condition", "coefficient_scale"
)]
```

The Stage 1 checkpoint additionally records whether execution used a supplied
BiocParallel backend, Pando's native map backend, or serial execution:

```r
step1$params$pando_parallel
```

## Compare the same reaction across conditions

After Stage 6, compare a fixed reaction/direction/medium target across
condition-by-broad-cell-type SuperCells:

```r
condition_stats <- rc_test_condition_reactions(
  result,
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  cell_types = "stem-cell_like",
  min_units = 5,
  p_adjust_method = "BH",
  p_adjust_scope = "celltype_contrast_medium"
)

condition_stats$omnibus
condition_stats$pairwise

reverse_report <- rc_report_condition_directions(
  result,
  target_directions = "reverse"
)

p <- rc_plot_condition_reaction(
  result,
  reaction_id = "MAR06231",
  cell_type = "epithelial_like",
  target_direction = "reverse",
  medium_scenario = "high_glucose",
  condition_col = "dataset",
  celltype_col = "epithelial_or_stem",
  conditions = c("control_24hr", "JQ1_24hr", "MS177_24hr"),
  annotation_p = "p_adj"
)
```

The omnibus test is Kruskal-Wallis; pairwise contrasts use Wilcoxon tests with
explicit multiplicity scope. The plot retains one point per metacell and
adjusted-P-value significance brackets. These are metacell-level descriptive
comparisons unless biological replicate-level inference is supplied separately.
See [Condition-associated reaction statistics](docs/condition-reaction-statistics.md).

## Restartable stages

- `rc_regcompass_step_grn()` — Pando condition GRNs and fit contracts.
- `rc_regcompass_step_metacells()` — condition-level multimodal metacells.
- `rc_regcompass_step_meta_modules()` — complete-GPR cores and biological catalogue.
- `rc_regcompass_step_layer1()` — RNA and RNA+ATAC reaction support.
- `rc_regcompass_step_layer2()` — shared structural model and directional LPs.
- `rc_regcompass_step_results()` — annotations, rankings, provenance, and contrasts.

## Tutorials

- [Level 1: minimal one-shot run](docs/tutorial-01-quick-start.md)
- [Level 2: stepwise run and audit](docs/tutorial-02-stepwise-audit.md)
- [Level 3: restart and sensitivity](docs/tutorial-03-advanced-restart.md)
- [Level 4: targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Level 5: condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Condition-associated reaction statistics](docs/condition-reaction-statistics.md)
- [Pando condition-comparable contract](docs/condition-comparable-grn.md)
- [Condition-comparability safeguards](docs/condition-comparability-safeguards.md)
- [Public API index](docs/functions.md)
