# Tutorial 2: restartable workflow

Each stage writes an RDS checkpoint to its output directory.

## Stage 1: regulatory evidence

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 6L, progressbar = TRUE)
}

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
      tf_cor = 0.05,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  parallel = TRUE,
  BPPARAM = upstream_bp
)
```

Stage 1 resolves the route independently for each retained cell type:

- at least two retained conditions: common-dictionary condition GRN;
- one retained condition: standard Pando;
- no retained stratum: excluded.

The same `pando_infer_args` list can be used for all routes. Condition-only arguments (`padj_threshold`, `rank_action`, `min_residual_df`, and layer controls) are disabled before standard Pando is called. Standard-model controls (`method`, `alpha`, `scale`, and related model arguments) are disabled for the fixed common-dictionary condition model. Unknown argument names still fail before model fitting.

Candidate discovery uses `abs(tf_target_cor) > 0.05` and `abs(peak_target_cor) > 0.05`. Final penalty entry requires `estimable == TRUE`, `padj < 0.05`, `abs(corr) >= 0.05`, and `abs(estimate) >= 0.05`.

When at least two cell-type jobs are available, Stage 1 distributes those jobs through `BPPARAM`. Every worker runs its own Pando job serially, so nested worker pools are not created. With one standard-Pando job, its existing target-level path may be used. A single condition-GRN job remains serial because pooled discovery, condition-specific discovery, dictionary freezing, and condition refits are treated as one coordinated contract.

```r
step1$grn_result$cell_type_analysis_mode
step1$grn_result$pando_execution_plan
step1$grn_result$pando_infer_argument_routing
```

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
    min_cells_per_stratum = 300L,
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

The default remains compact add-only FASTCORE completion:

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
      model_completion = "fastcore",
      completion_time_limit = 600,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 2000,
      strict = TRUE
    )
  )
)
```

To run the pinned Python CORDA2 semantics exactly for `met_prod = NULL`:

```r
layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 12L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 12L, progressbar = TRUE)
}

step5_corda2 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_corda2",
  model_mode = "meta_module_gem",
  parallel = TRUE,
  BPPARAM = layer2_bp,
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "corda2",
      completion_time_limit = 3000,
      strict = TRUE,
      corda2_redundancies = 3L,
      corda2_penalty_factor = 100,
      corda2_support = 5L,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

The reference is `resendislab/corda` commit
`c02e06d50606bf93f23d8f2e6d6ade0e996ca70e`. The implementation follows that
source as written rather than applying the previous PR-specific corrections.

CORDA2 maps merged core reactions to high confidence (`3`), non-core module reactions to medium confidence (`2`), optional high-evidence outside-module reactions to low confidence (`1`), low-evidence reactions to absent confidence (`-1`) and remaining reactions to unknown confidence (`0`). This mapping prepares the Python algorithm input; both forward and reverse variables initially receive the reaction confidence.

The exact source behavior is:

1. create both forward and reverse variables for every reaction, including closed directions;
2. normalize open reaction bounds to `UPPER = 1e6` using the active solver feasibility tolerance;
3. assess confidence-`3` targets serially and promote their associated variables;
4. assess remaining confidence-`1/2` targets serially, count associated absent variables, and promote those meeting `corda2_support`;
5. block remaining absent variables;
6. for each remaining confidence-`1/2` variable, minimize a positive coefficient and promote it only when the minimum objective is greater than `tflux = 1`;
7. block the remaining confidence-`1/2` variables, convert unknown variables to absent, and run one final dependency pass for the retained targets.

During target assessment, only the target variable receives `lb = max(1, lb)` and `ub = 1e6`. The opposite reversible variable is **not** closed. Costs for both directions of a reaction are assigned from the current forward-variable confidence, matching the Python source even after directional confidences diverge.

The only CORDA2 constructor controls exposed are:

- `corda2_redundancies` (`n`, default `3`);
- `corda2_penalty_factor` (`penalty_factor`, default `100`);
- `corda2_support` (`support`, default `5`).

`CI = 1.01`, `tflux = 1`, and `UPPER = 1e6` are fixed source constants. Do not supply `corda2_cost_increase`, `corda2_target_flux`, `corda2_flux_tolerance`, or `corda_seed`; these are rejected. `fastcore_epsilon` and `max_support_reactions` remain FASTCORE-only controls.

The parent model is the complete GEM after applying the requested medium bounds. CORDA2 does not run FASTCC, generic reaction-role pruning, or an added global feasibility precheck. Final models restore the original reaction bounds.

Only included core-reaction directions that reach the fixed `tflux = 1` on restored bounds enter downstream scoring. The Python build records impossible targets and completes, so `strict = TRUE` does not convert those results into a reconstruction error; they remain diagnostics and are excluded from the score matrices.

Each CORDA2 model is serial internally to preserve Python target order, confidence mutation order, and persistent solver state. Parallelism is used only across independent cell-type-by-medium model instances. A few models therefore will not use all workers during structural reconstruction. Model caches are written under `model_cache/meta_module_gem/corda2`.

See `docs/layer2-corda.md` for the full equations, line-by-line source correspondence, output contracts, audit fields, and Python-oracle CI design.

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

To restart, load the last valid checkpoint and rerun only later stages. Do not combine checkpoints created from different cell sets, GEMs, media, metadata columns, or Layer 2 completion methods.
