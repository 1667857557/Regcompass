# RegCompassR

RegCompassR 1.7.0 provides one canonical RNA+ATAC workflow:

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells with sample composition retained
→ one Pando GRN per condition × cell type
→ complete-GPR core reactions
→ core-reaction subsystem + KEGG/Reactome + master-Rhea expansion
→ local FASTCORE feasibility completion
→ one shared union-GEM and one shared extracellular medium
→ Pando-coefficient-weighted ATAC regulation integrated into RNA support
→ directional COMPASS-like minimum-penalty scoring
→ descriptive comparison between conditions
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

## Quick start

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
  metacell_args = list(
    gamma = 150,
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
      peak_cor = 0.05,
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
  )
)
```

The v1.7.0 canonical path requires `fragment_files = FALSE` and aggregates the
existing ATAC peak-count assay. Each biological sample must map to one condition.
With `strict_biological_defaults = TRUE`, each condition must contain at least
two biological samples. Cells are pooled by condition and cell type before
SuperCell2, but original sample membership is retained in
`result$metacells$sample_composition` and the corresponding output tables.

Condition-pooled metacells are descriptive pseudo-observations, not independent
biological replicates. The package does not perform biological-sample-level
significance testing on those metacells.

## Core multiome calculation

RNA is converted to zero-preserving bounded gene support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficients are learned from RNA+ATAC within each condition and cell
type. For per-metacell scoring, coefficient sign and magnitude weight robustly
standardized **peak accessibility only**. Metacell TF RNA is not multiplied into
the regulatory state, so target-gene RNA enters direct metacell support once.
The resulting modifier is bounded as \(R_{g,u}\in[-1,1]\).

The regulatory modifier changes RNA support on the support log-odds scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

This preserves zero RNA support, keeps values in `[0,1]`, increases support under
positive regulation and decreases it under negative regulation. Because the
Pando coefficients are fitted on the same pooled dataset, they are learned
parameters rather than independent validation evidence; external fitting or
cross-fitting is required for a fully independent regulatory layer.

Protein complexes use the normalized Boltzmann soft-min AND rule with
`tau = 0.20`:

\[
C_{complex}=-\tau\log\left(\frac{1}{n}\sum_{i=1}^{n}
\exp\left[-C_i/\tau\right]\right).
\]

Isozymes are added and no gene-promiscuity weighting is applied. Reaction
expression is converted to one COMPASS-like cost:

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

There is no independent Pando reaction-confidence penalty, Q95 calibration,
confidence-alignment matrix, or `penalty_weights` term in the canonical model.

## Structural model and LP

A reaction is core only when at least one complete GPR isozyme group is present.
Biological meta-module membership is then expanded only through:

1. the subsystem of each core reaction;
2. shared KEGG or Reactome reaction identifiers;
3. the same master Rhea identifier.

No reaction is added merely because it shares a metabolite with an included
reaction. There is no metabolite-neighbour or one-hop expansion API. Local
FASTCORE is the only stage that may add non-annotated reactions, and those
reactions are recorded separately as feasibility support rather than biological
meta-module members.

All conditions use the same union-GEM, stoichiometric matrix, bounds, medium,
target reactions and target-flux fraction. For each target direction, the solver
first obtains maximum feasible target flux, then constrains the target to at
least `omega × vmax` and minimizes the network-wide weighted absolute flux.

Primary outputs are:

- `result$metacells`: pooled metacells, membership and biological-sample composition;
- `result$layer1`: RNA support, ATAC-derived modifier, multiome gene support and `reaction_expression`;
- `result$grn_meta_modules`: annotation-defined biological membership, local FASTCORE support and shared union-GEM membership;
- `result$microcompass`: raw minimum penalties, feasibility and directional target diagnostics;
- `result$condition_summary` and `result$condition_contrast`: descriptive within-cell-type condition comparisons.

For the exact architecture and equations, see
[`docs/v1.7.0-condition-pooled-architecture.md`](docs/v1.7.0-condition-pooled-architecture.md).
