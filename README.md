# RegCompassR

RegCompassR 1.8.8 implements a shared-background regulatory–metabolic workflow for paired single-cell RNA+ATAC data.

## Canonical architecture

```text
all conditions within one cell type
→ one validated Pando structural TF–peak–metabolic-gene candidate universe
→ condition-balanced multitask elastic-net model
→ global GRN backbone + condition-specific deviations
→ condition-stratified bootstrap stability
→ condition-specific sub-GRNs and regulated metabolic genes
→ complete-GPR condition core reactions
→ ordered subsystem / KEGG–Reactome / master-Rhea expansion
→ one merged reaction catalogue
→ one medium-specific union GEM reused by every condition and metacell
→ RNA+ATAC reaction penalties
→ directional COMPASS-like LP scoring
```

For edge \(e=(TF,peak,target)\), the condition coefficient is

\[
\theta_{e,c}=\beta_e+\delta_{e,c},
\qquad
\sum_c\delta_{e,c}=0.
\]

`beta` is the cell-type-wide backbone coefficient. `delta` is the symmetric condition deviation. All conditions use the same edge dictionary, feature scaling, penalty structure, and stoichiometric model.

A condition-specific metabolic gene is supported when at least one bootstrap-stable active edge targets it. A reaction is a condition core only when one complete GPR branch is contained in that supported gene set. Positive and negative stable edges both establish that a gene is regulated.

## Installation

The validated default profile remains Seurat v4:

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@supercell-2.0"
)
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/Regcompass")
```

Pando 1.1.3 or later is required. RegCompass validates the version-2 Pando design fingerprint before fitting the multitask GRN.

SeuratObject/Seurat 5.x with Signac 1.12–1.x is also accepted. See [Seurat compatibility](docs/seurat-compatibility.md).

## Required metadata

The canonical workflow requires only:

```text
condition_col
celltype_col
```

No biological-sample column is accepted or interpreted by `rc_run_regcompass()`, `rc_regcompass_step_grn()`, or `rc_regcompass_step_metacells()`.

## Minimal complete run

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

  # Stage 1: shared candidate background + condition sub-GRNs
  grn_mode = "multitask_shared_backbone",
  pando_args = list(
    min_cells = 300,
    pando_design_args = list(
      peak_to_gene_method = "Signac",
      min_tf_detection = 0.01,
      min_peak_detection = 0.01,
      min_target_detection = 0.01
    )
  ),
  multitask_args = list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 2,
    lambda_rule = "lambda.1se",
    nfolds = 5,
    n_bootstrap = 100,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    seed = 12345L
  ),

  # Stage 2: condition-only, label-guided metacells
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30,
    seed = 12345L,
    min_cells_per_stratum = 300,
    min_metacell_size = 10,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),

  # Stage 4
  layer1_args = list(
    regulatory_alpha = 1,
    gpr_and_method = "min"
  ),

  # Stage 5
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
  upstream_workers = 6,
  layer2_workers = 30
)
```

The package default is `n_bootstrap = 50L`; `100` is shown for a final analysis with a more stable estimate of selection probabilities.

When `pfm` is omitted, RegCompass loads `data("motifs", package = "Pando")`. Unless `pando_args$pando_initiate_args$regions` is supplied, regulatory-region defaults are:

```text
human: phastConsElements20Mammals.UCSC.hg38 ∪ SCREEN.ccRE.UCSC.hg38
mouse: phastConsElements20Mammals.UCSC.hg38
```

## Core calculations

For a candidate edge \(e=(t,p,g)\), the Pando-style predictor is

\[
x_{e,u}=T_{t,u}A_{p,u}.
\]

Target expression and every TF-by-ATAC predictor are centred within condition. Edge scales are estimated once across all conditions of the cell type. Observation weights are proportional to \(1/n_c\), so each condition contributes the same total regression loss.

For bootstrap replicate \(b\), every condition is resampled with replacement at its original cell count. The resampled target and predictors are then re-centred within condition and fitted at the full-data selected \(\lambda\). This gives

\[
\Pi_{e,c}=\frac{1}{B}\sum_{b=1}^{B}
I\!\left(\theta^{(b)}_{e,c}\ne 0\right),
\]

and conditional sign stability

\[
\rho_{e,c}=
\left|
\frac{
\sum_b I(\theta^{(b)}_{e,c}\ne0)
\operatorname{sign}(\theta^{(b)}_{e,c})
}{
\sum_b I(\theta^{(b)}_{e,c}\ne0)
}
\right|.
\]

The reliability-weighted coefficient used by Layer 1 is

\[
\widetilde\theta_{e,c}
=
\widehat\theta_{e,c}\Pi_{e,c}\rho_{e,c}.
\]

Layer 1 converts these coefficients to an ATAC-only modifier. TFs sharing the same measured peak are signed-summed before the peak is projected. One target-specific denominator is shared across conditions, so weak condition networks are not independently rescaled to the same magnitude as strong networks. A missing active condition edge gives \(R=0\), and the log-odds update returns exact RNA-only support.

See [Shared-background multitask GRN mathematics and object contracts](docs/multitask-shared-grn.md).

## Inspectable stages

- `rc_regcompass_step_grn()`: build one validated Pando candidate background per cell type, fit global plus condition coefficients, and calculate bootstrap stability.
- `rc_regcompass_step_metacells()`: construct condition-only multimodal metacells with cell type used as the SuperCell2 label.
- `rc_regcompass_step_meta_modules()`: map condition sub-GRN targets to complete-GPR cores and biological reaction modules.
- `rc_regcompass_step_layer1()`: calculate RNA support, ATAC regulatory modifiers, and reaction expression.
- `rc_regcompass_step_layer2()`: build one shared medium-specific union GEM and run directional LP scoring.
- `rc_regcompass_step_results()`: assemble rankings, annotations, provenance, and condition contrasts.
- `rc_regcompass_step_target_union()`: score additional mapped targets in the exact cached union GEM.

## Main outputs

```r
result$grn$tf_peak_gene_candidates
result$grn$tf_peak_gene_global
result$grn$tf_peak_gene_condition_all
result$grn$tf_peak_gene_significant
result$grn$condition_target_genes
result$grn$stability_diagnostics
result$condition_grn_meta_modules$supported_metabolic_genes
result$condition_grn_meta_modules$core_gene_reaction
result$merged_grn_meta_modules$merged_core_reactions
result$merged_grn_meta_modules$merged_reaction_membership
result$reaction_ranking
result$condition_contrast
result$microcompass$model_cache_summary
```

## Legacy mode

Independent condition-by-cell-type Pando fitting remains available for reproducibility:

```r
result_legacy <- rc_run_regcompass_one_shot(
  ...,
  grn_mode = "legacy_condition_pando",
  pando_args = list(
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0.01,
      adjust_method = "fdr",
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1
  )
)
```

Legacy Pando thresholds belong only to that mode. They are not interpreted as significance thresholds in the multitask model.

## Tutorials

- [Quick start](docs/tutorial-01-quick-start.md)
- [Stepwise audit](docs/tutorial-02-stepwise-audit.md)
- [Advanced restart and diagnostics](docs/tutorial-03-advanced-restart.md)
- [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Condition differential analysis](docs/tutorial-05-condition-differential-analysis.md)
- [Stage contracts](docs/stage-interface-contracts.md)
- [Medium presets](docs/medium-presets.md)
