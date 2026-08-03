from pathlib import Path

ROOT = Path(".")

files = {
"README.md": r'''# RegCompassR

RegCompassR integrates paired single-cell RNA and ATAC regulatory evidence with cell-type-specific metabolic reaction analysis.

## Installation

The package uses the latest default branches of the companion repositories.

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

No companion commit or release is pinned in `DESCRIPTION`.

## Required input

The input is a paired-cell Seurat object containing:

- RNA and ATAC count assays for the same cells;
- a broad cell-type metadata column;
- RNA PCA and ATAC LSI reductions;
- genome-compatible peak coordinates;
- an optional condition column.

Stage 1 retains condition-by-cell-type strata with at least 300 paired cells. Within each cell type, two or more retained conditions use the common-dictionary condition GRN; one retained condition uses standard Pando automatically.

## One-shot workflow

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  )
)
```

## Documentation

- [Quick start](docs/tutorial-01-quick-start.md)
- [Restartable stages](docs/tutorial-02-stepwise-audit.md)
- [Targeted reaction remapping](docs/tutorial-04-targeted-reaction-remapping.md)
- [Condition-level results](docs/tutorial-05-condition-differential-analysis.md)
- [Mathematical specification](docs/mathematical-model.md)
- [Function index](docs/functions.md)
''',

"docs/tutorial-01-quick-start.md": r'''# Tutorial 1: quick start

## Input requirements

Use a paired-cell Seurat object with RNA and ATAC counts, broad cell-type metadata, RNA PCA, ATAC LSI, and genome-compatible peak coordinates.

Stage 1 uses a fixed minimum of 300 paired cells. For each cell type:

- at least two retained conditions: common-dictionary condition GRN;
- one retained condition: standard Pando;
- no retained condition stratum: excluded.

## Install current companion repositories

```r
remotes::install_github("1667857557/Pando_regcompass")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
remotes::install_github("1667857557/Regcompass")
```

## Prepare GEM and medium

```r
library(RegCompassR)

gem <- rc_prepare_gem(
  species = "human",
  version = "2.0.0",
  source = "bundled"
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Supported built-in scenarios include:

| Scenario | Use |
|---|---|
| `normal_human_plasma` | Human plasma-like background |
| `mouse_plasma` | Mouse plasma/interstitial-fluid background |
| `high_glucose` | Plasma-like background with increased glucose |
| `low_glucose` | Plasma-like background with reduced glucose |
| `high_lactate` | Plasma-like background with increased lactate |
| `low_lactate` | Plasma-like background with reduced lactate |
| `low_glutamine` | Plasma-like background with reduced glutamine |
| `custom` | User-supplied exchange bounds or metabolite availability |

Custom reaction bounds:

```r
custom_medium <- data.frame(
  medium_scenario_id = "measured_medium",
  exchange_reaction_id = c("EX_glc_D_e", "EX_gln_L_e"),
  lb = c(-0.20, -0.10),
  ub = c(1, 1),
  available = TRUE
)

medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "custom",
  species = "human",
  custom_medium = custom_medium
)
```

## Run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  gem = gem,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  fragment_files = FALSE,
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    k.knn = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L,
    overwrite = FALSE
  ),
  layer1_args = list(
    gpr_and_method = "min",
    gene_half_saturation = 1
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

## Main outputs

```r
result$grn$cell_type_analysis_mode
result$grn$condition_fit_status
result$layer1$gene_regulatory_modifier
result$microcompass$penalty
result$reaction_ranking
result$condition_contrast
```

The mathematical definitions are maintained only in [mathematical-model.md](mathematical-model.md).
''',

"docs/tutorial-02-stepwise-audit.md": r'''# Tutorial 2: restartable workflow

Each stage writes an RDS checkpoint to its output directory.

## Stage 1: regulatory evidence

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  )
)
```

`step1$grn_result$cell_type_analysis_mode` reports the route used by each cell type.

## Stage 2: metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
  metacell_args = list(
    rna_reduction = "pca",
    rna_dims = 1:30,
    atac_reduction = "lsi",
    atac_dims = 2:30,
    gamma = 30L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  )
)
```

Stage 2 reproduces the exact Stage 1 cell set. One WNN graph is built per broad cell type; final metacells remain condition-pure.

## Stage 3: biological reaction catalogue

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)
```

## Stage 4: Layer 1 reaction support

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  gpr_and_method = "min",
  gene_half_saturation = 1
)
```

## Stage 5: Layer 2 metabolic scoring

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
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
  )
)
```

## Stage 6: result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "run/06_results",
  species = "human"
)
```

To restart, load the last valid checkpoint and rerun only later stages. Do not combine checkpoints created from different cell sets, GEMs, media, or metadata columns.
''',

"docs/tutorial-04-targeted-reaction-remapping.md": r'''# Tutorial 4: targeted reaction remapping

Use targeted remapping after a completed RegCompass run when a focused reaction list is required.

## Define targets

```r
targets <- data.frame(
  reaction_id = c("MAR00123", "MAR00456"),
  target_direction = c("forward", "both"),
  stringsAsFactors = FALSE
)
```

Reaction identifiers must map to the selected GEM. Keep direction explicit when reversible reactions are interpreted separately.

## Run targeted remapping

```r
targeted <- rc_regcompass_targeted_reactions(
  result = result,
  gem = gem,
  target_reactions = targets,
  outdir = "RegCompass_targeted",
  solver = "highs"
)
```

The function reuses compatible cached structural models and medium settings. A different GEM, medium table, cell-type catalogue, direction setting, or cache checksum requires rebuilding the affected model.

## Inspect outputs

```r
targeted$target_reactions
targeted$reaction_ranking
targeted$condition_summary
targeted$condition_contrast
targeted$model_cache_summary
```

Targeted remapping does not refit Pando or rebuild metacells.
''',

"docs/tutorial-05-condition-differential-analysis.md": r'''# Tutorial 5: condition-level reaction results

Condition comparisons use metacells within the same cell type, reaction, direction, and medium.

## Inspect available groups

```r
unique(result$reaction_comparison_by_metacell[, c(
  "cell_type", "condition", "medium", "direction"
)])
```

Cell types routed to standard Pando because only one condition was retained contribute reaction rankings but cannot produce a within-cell-type condition contrast.

## Reaction ranking

```r
ranking <- result$reaction_ranking
ranking <- ranking[
  ranking$cell_type == "T_cell" &
  ranking$medium == "normal_human_plasma",
]
```

## Pairwise condition contrast

```r
contrast <- result$condition_contrast
contrast <- contrast[
  contrast$cell_type == "T_cell" &
  contrast$condition_1 == "Control" &
  contrast$condition_2 == "Treatment",
]
```

Review the number of metacells per condition and the feasible/evaluated flags before interpreting a comparison.

## RNA-only control

```r
result$rna_only_control_summary
result$rna_only_control_contrast
```

The RNA-only route uses the same structural models and media and is intended as an interpretation control.

## Export

```r
write.table(
  ranking,
  file = "reaction_ranking.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

Reported statistical tests describe within-dataset metacell separation. Donor-level inference requires donor-aware biological replication outside this workflow.
''',

"vignettes/regcompass-workflow.Rmd": r'''---
title: "RegCompass workflow"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{RegCompass workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE, collapse = TRUE, comment = "#>")
```

# Input

Use a paired RNA+ATAC Seurat object with broad cell-type metadata, RNA PCA, ATAC LSI, and an optional condition column.

Stage 1 routes each cell type independently:

- two or more retained conditions: condition GRN;
- one retained condition: standard Pando.

# One-shot run

```r
result <- rc_run_regcompass_one_shot(
  object = A,
  outdir = "RegCompass_result",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  species = "human",
  condition_col = "condition",
  celltype_col = "cell_type",
  pando_args = list(
    min_cells = 300L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  metacell_args = list(
    rna_reduction = "pca",
    atac_reduction = "lsi",
    gamma = 30L
  )
)
```

# Stepwise run

```r
step1 <- rc_regcompass_step_grn(...)
step2 <- rc_regcompass_step_metacells(..., grn = step1)
step3 <- rc_regcompass_step_meta_modules(
  grn = step1, metacells = step2, gem = gem, ...
)
step4 <- rc_regcompass_step_layer1(
  grn = step1, metacells = step2, meta_modules = step3, gem = gem, ...
)
step5 <- rc_regcompass_step_layer2(
  layer1 = step4, meta_modules = step3, gem = gem, ...
)
result <- rc_regcompass_step_results(
  grn = step1, metacells = step2, meta_modules = step3,
  layer1 = step4, layer2 = step5, gem = gem, ...
)
```

# Outputs

```r
result$grn$cell_type_analysis_mode
result$reaction_ranking
result$condition_contrast
result$rna_only_control_contrast
```

The mathematical specification is maintained only in `docs/mathematical-model.md`.
''',

"docs/functions.md": r'''# Function index

## Complete workflow

- `rc_run_regcompass_one_shot()`: prepare defaults and run the complete workflow.
- `rc_run_regcompass()`: run the complete workflow with explicit stage arguments.

## Restartable stages

- `rc_regcompass_step_grn()`: filter cells and fit cell-type-specific Pando routes.
- `rc_regcompass_step_metacells()`: construct condition-pure SuperCell metacells.
- `rc_regcompass_step_meta_modules()`: build condition-by-cell-type biological reaction catalogues.
- `rc_regcompass_step_layer1()`: combine RNA and regulatory support and apply GPR rules.
- `rc_regcompass_step_layer2()`: build structural models and score directional reactions.
- `rc_regcompass_step_results()`: assemble annotations, rankings, and contrasts.

## GEM and media

- `rc_prepare_gem()`: load and prepare a supported GEM.
- `rc_validate_gem()`: validate a prepared GEM.
- `rc_make_medium_scenarios()`: construct built-in or custom medium tables.

## Results

- `rc_regcompass_targeted_reactions()`: rerun scoring for a focused reaction list.
- `rc_regcompass_condition_contrast()`: extract condition-level comparisons where available.
- `rc_export_microcompass()`: export Layer 2 matrices and diagnostics.

Use the generated Rd help for complete argument definitions. Mathematical definitions are in [mathematical-model.md](mathematical-model.md).
''',

"NEWS.md": r'''# RegCompassR 2.4.0

- Uses the latest default branches of Pando_regcompass and SuperCell_Seurat_V4 without fixed revisions.
- Routes each retained cell type independently to condition GRN or standard Pando according to its retained condition count.
- Uses the current Pando condition-fit API.
- Removes retired projection and penalty fields from Layer 1, Layer 2, documentation, and result schemas.
- Consolidates mathematical definitions in `docs/mathematical-model.md`.
- Reduces user tutorials and continuous integration to the current supported workflow.
'''
}

for name, content in files.items():
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.rstrip() + "\n")

remove_docs = [
    "docs/tutorial-03-mathematical-model.md",
    "docs/condition-comparable-grn.md",
    "docs/condition-comparability-safeguards.md",
    "docs/canonical-defaults.md",
    "docs/workflow.md",
    "docs/stage-interface-contracts.md",
    "docs/target-union-scoring.md",
    "docs/run-modes-and-stepwise-workflow.md"
]
for name in remove_docs:
    path = ROOT / name
    if path.exists():
        path.unlink()

workflow = ROOT / ".github/workflows/fastcore-checks.yaml"
workflow.write_text(r'''name: RegCompassR package check

on:
  pull_request:
  push:
    branches: [Main]

permissions:
  contents: read

jobs:
  check:
    runs-on: ubuntu-22.04
    timeout-minutes: 90
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          r-version: '4.5'
          use-public-rspm: true
      - uses: r-lib/actions/setup-pandoc@v2
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          dependencies: '"hard"'
          extra-packages: |
            any::rcmdcheck
            github::1667857557/Pando_regcompass
            github::1667857557/SuperCell_Seurat_V4
      - name: Parse R sources
        shell: Rscript {0}
        run: |
          files <- c(
            list.files("R", pattern = "[.]R$", full.names = TRUE),
            list.files("tests/testthat", pattern = "[.]R$", full.names = TRUE)
          )
          for (file in files) parse(file = file)
      - name: Check package
        shell: Rscript {0}
        env:
          _R_CHECK_FORCE_SUGGESTS_: 'false'
        run: rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
''')

for path in (ROOT / ".github/workflows").glob("*"):
    if path.name != "fastcore-checks.yaml":
        path.unlink()

for name in [
    "tests/testthat/test-current-condition-contract.R",
    "tests/testthat/test-pando-common-dictionary-runtime.R"
]:
    path = ROOT / name
    if path.exists():
        path.unlink()

strict = ROOT / "tests/testthat/test-strict-common-dictionary-penalty.R"
if strict.exists():
    text = strict.read_text()
    marker = 'test_that("package metadata pins the audited Pando head"'
    start = text.find(marker)
    if start >= 0:
        strict.write_text(text[:start].rstrip() + "\n")

replacements = {
    "reaction_expression_condition_full_oof": "reaction_expression",
    "reaction_expression_common_oof": "reaction_expression",
    "gene_projection_condition_full_oof": "gene_projection",
    "gene_projection_common_oof": "gene_projection",
    "gene_support_condition_full_oof": "gene_support_multiome",
    "gene_support_common_oof": "gene_support_multiome",
    "gene_regulatory_modifier_condition_full_oof": "gene_regulatory_modifier",
    "gene_regulatory_modifier_common_oof": "gene_regulatory_modifier",
    "penalty_condition_full_oof": "penalty",
    "penalty_common_oof": "penalty",
    "score_condition_full_oof_display_only": "score",
    "score_common_oof_display_only": "score",
    "score_rna_only_display_only": "score_rna_only",
    "condition_full_oof": "primary",
    "common_oof": "primary",
    "condition_unique_oof": "removed",
    "condition_unique_increment": "removed_increment",
    "min_abs_estimate": "absolute_estimate_threshold",
    "min_model_rsq": "model_quality_threshold"
}
for root in [ROOT / "docs", ROOT / "vignettes", ROOT / "man", ROOT / "tests"]:
    if not root.exists():
        continue
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        for old, new in replacements.items():
            text = text.replace(old, new)
        path.write_text(text)

for path in list((ROOT / "docs").rglob("*.md")) + \
            list((ROOT / "vignettes").rglob("*.Rmd")) + [ROOT / "README.md"]:
    if path == ROOT / "docs/mathematical-model.md":
        continue
    text = path.read_text()
    for token in ("\\[", "\\]", "\\(", "\\)", "$$"):
        if token in text:
            raise RuntimeError(
                f"mathematical markup remains outside dedicated file: {path}"
            )

for token in (
    "condition_full_oof", "common_oof", "condition_unique_oof",
    "min_abs_estimate", "min_model_rsq"
):
    hits = []
    for root in [ROOT / "R", ROOT / "docs", ROOT / "vignettes",
                 ROOT / "man", ROOT / "tests", ROOT / "README.md"]:
        paths = [root] if root.is_file() else list(root.rglob("*"))
        for path in paths:
            if not path.is_file():
                continue
            try:
                text = path.read_text()
            except UnicodeDecodeError:
                continue
            if token in text:
                hits.append(str(path))
    if hits:
        raise RuntimeError(f"retired token {token} remains in: {hits}")

scripts = ROOT / ".github/scripts"
if scripts.exists():
    for path in scripts.glob("*"):
        path.unlink()
