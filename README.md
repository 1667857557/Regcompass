# RegCompassR

RegCompassR 1.7.0 provides one condition-pooled RNA+ATAC workflow with two execution modes: a one-shot convenience runner and inspectable, restartable public stages.

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells with sample composition retained
→ zero-count ATAC peaks excluded
→ one Pando GRN per condition × cell type
→ complete-GPR core reactions
→ core-reaction subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ shared union meta-module GEM or shared full GEM
→ Pando-coefficient-weighted ATAC regulation integrated into RNA support
→ directional COMPASS-like minimum-penalty scoring
→ within-condition reaction ranking
→ optional descriptive pairwise comparison when multiple conditions are present
```

## Installation

```r
install.packages("remotes")
remotes::install_version("SeuratObject", "4.1.4", upgrade = "never")
remotes::install_version("Seurat", "4.4.0", upgrade = "never")
remotes::install_version("Signac", "1.11.0", upgrade = "never")
remotes::install_github(
  "1667857557/SuperCell_Seurat_V4@supercell-2.0",
  upgrade = "never"
)
remotes::install_github("1667857557/Pando_regcompass", upgrade = "never")
remotes::install_github("1667857557/Regcompass", upgrade = "never")
```

RegCompassR validates this exact Seurat stack when the package is loaded.

## Execution mode A: one-shot

Use one-shot mode after the analysis parameters have been established.

```r
library(RegCompassR)
library(Pando)
library(BSgenome.Hsapiens.UCSC.hg38)

data(motifs, package = "Pando")
data(SCREEN.ccRE.UCSC.hg38, package = "Pando")

result <- rc_run_regcompass_one_shot(
  object = object,
  outdir = "RegCompass_result",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  fragment_files = FALSE,
  species = "human",
  gem_version = "2.0.0",
  medium_scenario = "normal_human_plasma",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  model_mode = "meta_module_gem",
  metacell_args = list(
    gamma = 20,
    min_cells_per_stratum = 500,
    min_metacell_size = 10
  ),
  pando_args = list(
    min_metacells = 10,
    pando_initiate_args = list(
      regions = SCREEN.ccRE.UCSC.hg38
    ),
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    )
  ),
  layer1_args = list(
    regulatory_alpha = 1,
    tau = 0.20,
    local_fastcore = TRUE
  ),
  layer2_args = list(
    target_direction = "both",
    solver = "highs"
  ),
  upstream_workers = 8,
  layer2_workers = 8
)
```

## Execution mode B: inspectable stepwise workflow

Use stepwise mode while selecting parameters, auditing intermediate biology or restarting an expensive analysis from a saved checkpoint.

```r
gem <- rc_prepare_gem(species = "human", version = "2.0.0")
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)

bp <- BiocParallel::SnowParam(workers = 8, type = "SOCK")

step1 <- rc_regcompass_step_metacells(
  object = object,
  outdir = "RegCompass_steps/01_metacells",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "cell_type",
  fragment_files = FALSE,
  metacell_args = list(
    gamma = 20,
    min_cells_per_stratum = 500,
    min_metacell_size = 10,
    BPPARAM = bp
  )
)

step2 <- rc_regcompass_step_meta_modules(
  metacells = step1,
  gem = gem,
  outdir = "RegCompass_steps/02_meta_modules",
  pfm = motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pando_args = list(
    min_metacells = 10,
    pando_initiate_args = list(regions = SCREEN.ccRE.UCSC.hg38),
    pando_infer_args = list(
      method = "glm",
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "fdr"
    ),
    padj_threshold = 0.05,
    min_model_rsq = 0.1
  ),
  layer1_args = list(local_fastcore = TRUE),
  parallel = TRUE,
  BPPARAM = bp
)

step3 <- rc_regcompass_step_layer1(
  metacells = step1,
  meta_modules = step2,
  gem = gem,
  outdir = "RegCompass_steps/03_layer1",
  regulatory_alpha = 1,
  tau = 0.20,
  gene_half_saturation = 1,
  parallel = TRUE,
  BPPARAM = bp
)

step4 <- rc_regcompass_step_layer2(
  layer1 = step3,
  meta_modules = step2,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "RegCompass_steps/04_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    omega = 0.95,
    solver = "highs"
  ),
  parallel = TRUE,
  BPPARAM = bp
)

result <- rc_regcompass_step_results(
  metacells = step1,
  meta_modules = step2,
  layer1 = step3,
  layer2 = step4,
  gem = gem,
  outdir = "RegCompass_steps/final",
  species = "human"
)
```

`parallel = FALSE` forces sequential execution in Steps 2–4. `BPPARAM = NULL` uses the default RegCompass backend, `BPPARAM = FALSE` forces base `lapply()`, and a `BiocParallelParam` object selects an explicit backend. Logical `TRUE` is not a valid `BPPARAM` value.

Each stage writes a complete RDS checkpoint and returns the object consumed by the next stage. Inspect `step1$pooled`, `step2$condition_modules`, `step2$global_modules`, `step3$gene_regulatory_modifier`, `step3$reaction_expression`, `step4$lp_diagnostics`, and the final reaction rankings before proceeding.

The detailed tutorial lists inspection commands, expected ranges, failure checks and the exact downstream stages that must be rerun after each parameter change:
[`docs/run-modes-and-stepwise-workflow.md`](docs/run-modes-and-stepwise-workflow.md).

The canonical path requires `fragment_files = FALSE` and aggregates the existing ATAC peak-count assay. Peaks with zero total counts are removed before shared TF-IDF normalization, and peaks with zero counts within a Pando condition-by-cell-type group are removed again before motif and GRN inference. The numbers retained and excluded are recorded in the ATAC normalization metadata and Pando status table.

Each biological sample must map to exactly one condition. One condition is sufficient for reaction ranking; two or more conditions additionally produce descriptive pairwise priority comparisons. Cells are pooled by condition and cell type before SuperCell2, while original sample membership remains available in `result$metacells$sample_composition`.

Condition-pooled metacells are descriptive pseudo-observations, not independent biological replicates. The package warns when a condition has fewer than two biological samples and does not perform biological-sample-level significance testing on metacells.

## Structural model modes

`model_mode = "meta_module_gem"` scores the shared union of all condition-specific biological meta-modules after local FASTCORE feasibility completion.

`model_mode = "full_gem"` scores the same targets and multiome penalties in the shared full GEM. Both modes use the same medium, bounds, target-flux fraction, Layer 1 evidence model, and descriptive ranking outputs.

## Core multiome calculation

RNA is converted to zero-preserving bounded gene support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficients are learned from RNA+ATAC within each condition and cell type. For per-metacell scoring, coefficient sign and magnitude weight robustly standardized **peak accessibility only**. Metacell TF RNA is not multiplied into the regulatory state, so target-gene RNA enters direct metacell support once. The resulting modifier is bounded as \(R_{g,u}\in[-1,1]\).

The regulatory modifier changes RNA support on the support log-odds scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

This preserves zero RNA support, keeps values in `[0,1]`, increases support under positive regulation and decreases it under negative regulation. Pando coefficients fitted on the same pooled dataset remain learned parameters rather than independent validation evidence.

Protein complexes use the normalized Boltzmann soft-min AND rule with `tau = 0.20`:

\[
C_{complex}=-\tau\log\left(\frac{1}{n}\sum_{i=1}^{n}
\exp\left[-C_i/\tau\right]\right).
\]

Isozymes are added and no gene-promiscuity weighting is applied. Reaction expression is converted to one COMPASS-like cost:

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

For expression-linked reactions, unmeasured reaction expression is set to zero before applying this formula. Therefore unmeasured expression and explicit zero expression both receive the strictest expression-linked LP penalty, while missingness remains available as a diagnostic flag.

Only structural exchange, demand, sink, and artificial-support reactions receive fixed structural costs. Transport and cofactor reactions with expression support remain governed by the same multiome reaction-expression cost. There is no independent Pando reaction-confidence penalty, Q95 calibration, confidence matrix, or `penalty_weights` term.

## Structural model, LP and ranking

A reaction is core only when at least one complete GPR isozyme group is present. Biological meta-module membership is expanded only through:

1. the subsystem of each core reaction;
2. shared KEGG or Reactome reaction identifiers;
3. the same master Rhea identifier.

No reaction is added merely because it shares a metabolite with an included reaction. Local FASTCORE is the only stage that may add non-annotated reactions, and those reactions are recorded separately as feasibility support.

All scored units use the same structural GEM for the selected mode, stoichiometric matrix, bounds, medium, target reactions and target-flux fraction. For each target direction, the solver first obtains maximum feasible target flux, then constrains the target to at least `omega × vmax` and minimizes the network-wide weighted absolute flux.

Raw LP penalties are used for comparisons of the same target between conditions. Cross-reaction priority uses

\[
\widetilde P_{r,u}=\frac{P^*_{r,u}}{\omega V^{max}_{r}},
\]

implemented as `penalty / (omega * vmax)`. This is the minimum evidence cost per unit required near-maximal target flux in the selected structural model. Reactions are ranked within each condition and cell type by ascending median normalized cost. This is a model-based priority, not a physical flux estimate.

Primary outputs are:

- `result$metacells`: pooled metacells, membership and sample composition;
- `result$layer1`: RNA support, ATAC-derived modifier, multiome gene support and reaction expression;
- `result$grn_meta_modules`: biological membership, local FASTCORE support and shared union membership;
- `result$microcompass`: raw minimum penalties, feasibility, `vmax` and directional target diagnostics;
- `result$reaction_ranking`: raw and flux-normalized reaction priorities within every condition and cell type;
- `result$condition_contrast`: descriptive pairwise comparisons when at least two conditions are present.

For the exact architecture and equations, see [`docs/v1.7.0-condition-pooled-architecture.md`](docs/v1.7.0-condition-pooled-architecture.md).
