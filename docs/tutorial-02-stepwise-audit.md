# Tutorial 2: restartable workflow

Each stage writes an RDS checkpoint to its output directory. Use the same paired-cell object, GEM, metadata columns and medium definitions when restarting later stages.

## 1. Parallel backend

```r
library(BiocParallel)

upstream_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 6L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 6L, progressbar = TRUE)
}

layer2_bp <- if (.Platform$OS.type == "windows") {
  SnowParam(workers = 30L, type = "SOCK", progressbar = TRUE)
} else {
  MulticoreParam(workers = 30L, progressbar = TRUE)
}
```

## 2. Regulatory evidence

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

A cell type with at least two retained conditions uses the common-dictionary condition GRN; this tutorial passes `tf_cor = 0.1` to that route. A cell type with one retained condition uses standard Pando. For standard Pando, the requested `tf_cor` is an effect-size floor and RegCompass computes an additional sample-size floor from the exact two-sided Pearson correlation test,

\[
r_{crit}(n) = \frac{t_{1-\alpha/2,n-2}}{\sqrt{t_{1-\alpha/2,n-2}^2+n-2}},\qquad \alpha=0.05,
\]

then passes `max(tf_cor, r_crit(n))` to `Pando::infer_grn()` for that cell type. Thus small standard-Pando groups require a stronger TF-target correlation solely because their null correlation distribution is wider. The condition-GRN route is not altered by this adaptive standard-Pando gate.

## 3. Metacells

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

One WNN graph is constructed per broad cell type. Final metacells remain condition-pure.

## 4. Reaction catalogue and Layer 1 support

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules"
)

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

## 5. Medium scenarios

Use `rc_make_medium_scenarios()` so every condition within one scenario receives identical exchange bounds.

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human"
)
```

Supported biological presets are:

```text
normal_human_plasma
mouse_plasma
high_glucose
low_glucose
high_lactate
low_lactate
low_glutamine
custom
```

Multiple scenarios can be supplied together:

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = c(
    "normal_human_plasma",
    "high_glucose",
    "low_glucose",
    "high_lactate",
    "low_lactate",
    "low_glutamine"
  ),
  species = "human"
)
```

For exact reaction-level custom bounds use `scenario = "custom"` and `custom_medium`. For metabolite-level composition use `scenario = NULL` and `custom_metabolites`. Provenance fields such as `background_reference_doi`, `background_validation_reference_doi` and `challenge_reference_doi` are retained. See `docs/medium-presets.md`.

## 6. Layer 2: default CORDA2 reconstruction

With `model_mode = "meta_module_gem"`, omitting `model_completion` now selects original MATLAB CORDA2 semantics.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  parallel = TRUE,
  BPPARAM = layer2_bp,
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      strict = TRUE,
      corda2_args = list(
        MCxNCthresh = 2,
        constraint = 1,
        constrainby = "val",
        om = 1e4,
        ci = 0.01
      ),
      corda_medium_confidence_threshold = 0.75,
      corda_negative_confidence_threshold = 0.10,
      corda_regulatory_weight = 0.20,
      corda_include_evidence_outside_modules = TRUE,
      corda_max_medium_confidence_reactions = Inf
    )
  )
)
```

The complete medium-constrained parent GEM is passed directly to CORDA2 without FASTCC pre-pruning. Final retained reactions recover the parent bounds for selected directions, including positive parent lower bounds. Layer 2 reaction costs use the COMPASS scale: missing expression and structural roles receive the maximum cost `1`.

CORDA2 reconstruction has no structural time limit. `model_params$completion_time_limit` is rejected on the default CORDA2 route so a long Human-GEM reconstruction cannot be silently truncated. The parameter remains available for supplementary non-CORDA2 completion such as FASTCORE.

With `parallel = TRUE`, `BPPARAM` supplies the Layer-2 worker budget and backend. CORDA2 does not parallelize different mathematical steps simultaneously. Instead, Step 1, Step 2.1, Step 2.2 and Step 3 remain strict barriers; directional targets inside the current step are distributed across the available workers, restored to original directional order, and only then are HC/MC/NC/OT states updated. Each step uses a fresh worker pool with one HiGHS thread per worker. The pool and worker-local solver engines are released and full garbage collection is requested before the next step starts. Step entry and completion are printed, and the step progress display reports `completed/total` plus `remaining` directional targets. The same target-level status is written to the Layer-2 task progress files with `scope = "corda2_stage"`.

Main controls:

- `target_direction`: `"both"`, `"forward"` or `"reverse"`;
- `solver`: `"highs"`, `"gurobi"` or `"glpk"`;
- `strict`: fail when required targets are not retained;
- `corda2_args`: original CORDA2 controls `MCxNCthresh`, `constraint`, `constrainby`, `om` and `ci`;
- `corda_*`: RegCompass evidence-to-confidence mapping controls;
- `completion_time_limit`: non-CORDA2 structural completion control; do not use it with CORDA2.

## 7. Supplementary structural modes

### FASTCORE

Set the completion method explicitly:

```r
step5_fastcore <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_fastcore",
  model_mode = "meta_module_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    model_params = list(
      model_completion = "fastcore",
      completion_time_limit = 1200,
      fastcore_epsilon = 1e-4,
      max_support_reactions = 3000,
      strict = TRUE
    )
  )
)
```

### Full GEM

Full-GEM mode skips context-specific reconstruction and applies only the medium exchange bounds before COMPASS-style directional scoring.

```r
step5_full <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2_full_gem",
  model_mode = "full_gem",
  layer2_args = list(
    target_direction = "both",
    solver = "highs",
    flux_threshold = 1e-8
  )
)
```

Do not supply CORDA2 or FASTCORE controls in `full_gem` mode.

## 8. Results

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

Useful Layer 2 provenance fields are `step5$params`, `step5$completion_contract`, `step5$model_cache_summary`, `step5$vmax_cache_diagnostics` and `step5$lp_diagnostics`.
