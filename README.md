# RegCompassR

RegCompassR 1.7.0 provides one canonical RNA+ATAC workflow:

```text
condition × cell type cells pooled across biological samples
→ SuperCell2 condition-level metacells
→ condition × cell type Pando GRNs
→ condition-specific GRN meta-modules and local FASTCORE completion
→ one shared union-GEM and one shared extracellular medium
→ TF–ATAC regulation integrated into gene support before GPR aggregation
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

The canonical v1.7.0 workflow requires `fragment_files = FALSE` and aggregates
the existing ATAC peak-count assay. Cells from all biological samples in the
same condition and cell type are supplied together to SuperCell2. Consequently,
a pooled metacell does not retain one biological-sample identity and the
canonical inference unit is the metacell.

## Core multiome calculation

RNA is converted to zero-preserving bounded gene support:

\[
C^{RNA}_{g,u}=\frac{x_{g,u}}{x_{g,u}+h}.
\]

Pando coefficients determine the sign and relative weight of TF–peak–gene
regulation. Peak accessibility and TF abundance are multiplied, centered and
scaled across all conditions within the same cell type, and aggregated to a bounded
regulatory modifier \(R_{g,u}\in[-1,1]\).

The regulatory modifier changes RNA support on the support log-odds scale:

\[
C^{MO}_{g,u}=
\frac{C^{RNA}_{g,u}2^{\alpha R_{g,u}}}
{1-C^{RNA}_{g,u}+C^{RNA}_{g,u}2^{\alpha R_{g,u}}}.
\]

This preserves zero RNA support, keeps all values in `[0,1]`, increases support
under positive regulation and decreases it under negative regulation.

Protein complexes use a Boltzmann minimum-biased AND rule with `tau = 0.20`.
Isozymes are added. No gene-promiscuity weighting is applied. The resulting
multiome reaction expression is converted to one COMPASS-like cost:

\[
p_{r,u}=\frac{1}{1+\log_2(1+E^{MO}_{r,u})}.
\]

Pando is therefore not added as an independent reaction-level penalty.

## Structural model and comparison

Each condition-specific GRN meta-module is completed locally with FASTCORE.
Completed modules from all conditions are deduplicated into one union-GEM. All
conditions use the same stoichiometric matrix, bounds, extracellular medium,
target reactions and target-flux fraction.

Primary outputs are:

- `result$layer1`: RNA support, TF–ATAC modifier, multiome gene support and reaction expression;
- `result$grn_meta_modules`: condition-specific modules and the shared union-GEM membership;
- `result$microcompass`: raw minimum penalties and directional target diagnostics;
- `result$condition_summary`: per-cell-type, per-condition median reaction penalties and support scores;
- `result$condition_contrast`: two-condition relative support differences within each cell type.

For the exact architecture and equations, see
[`docs/v1.7.0-condition-pooled-architecture.md`](docs/v1.7.0-condition-pooled-architecture.md).
