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

To run the original three-stage CORDA dependency-assessment algorithm:

```r
layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 12L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 12L, progressbar = TRUE)
}

step5_corda <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_corda",
  model_mode = "meta_module_gem",
  parallel = TRUE,
  BPPARAM = layer2_bp,
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "corda",
      completion_time_limit = 3000,
      strict = TRUE,
      corda_gamma = 1e5,
      corda_kappa = 1e-2,
      corda_epsilon = 1,
      corda_n = 5L,
      corda_p = 2L,
      corda_seed = 1L,
      corda_flux_tolerance = 1e-8,
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

CORDA maps merged core reactions to HC, non-core module reactions to flexible MC,
low-evidence reactions to NC and remaining reactions to OT. It then performs the
published HC-association, shared-NC/MC-feasibility and RE/OT-association stages.
The parent model is the complete GEM after applying the requested medium bounds;
CORDA does not apply FASTCC pruning or generic demand/sink/artificial-support
blocking. `fastcore_epsilon` and `max_support_reactions` are therefore FASTCORE
controls and should not be supplied to the CORDA route.

For a reversible reaction being assessed in one direction, the opposite split
direction is fixed to zero. This is equivalent to the original implementation,
which leaves the target reaction unsplit and constrains it to `+epsilon` or
`-epsilon`, and prevents a zero-net-flux forward/reverse self-cycle.

Only retained HC core directions that can reach `corda_epsilon` remain downstream
scoring targets, so the score and penalty matrix dimensions remain compatible
with the FASTCORE route. The optional temporary metabolite-task reactions from
the original MATLAB workflow are not exposed by the current Layer-2 input
contract.

When `BPPARAM` is omitted, the CORDA route creates one package-managed worker
pool and keeps it alive across all stages. With few cell-type-by-medium models,
target-direction-by-repeat tasks use the full pool; with many models, models run
in parallel. HiGHS uses a persistent native C++ solver and reuses the simplex
basis within each task block. Every worker restricts its internal LP solver to
one thread to avoid oversubscription.

See `docs/layer2-corda.md` for the exact confidence mapping, mathematical
transformation, original-code equivalence, acceleration strategy and audit
fields.

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
