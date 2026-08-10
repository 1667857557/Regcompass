# Tutorial 2: restartable workflow

Each stage writes an RDS checkpoint to its output directory. Use the same paired-cell object, GEM, metadata columns and medium definitions when restarting later stages.

## 1. One parallel worker cap

RegCompass selects the parallel backend automatically:

- Windows: `BiocParallel::SnowParam(type = "SOCK")`;
- Linux/macOS: `BiocParallel::MulticoreParam`.

Users set only one worker cap. The default is `10L`; increase it explicitly when more parallel capacity is desired:

```r
workers <- 60L
```

The effective worker count is always bounded by

```text
min(independent tasks, workers, max(1, detected logical CPUs - 2))
```

so RegCompass reserves two logical CPUs for the operating system/R controller. Condition-GRN candidate and fit phases can expose condition × cell-type jobs; standard Pando exposes broad-cell-type jobs; CORDA2 and directional LP phases can expose thousands of targets. Every dispatch uses only the workers it can actually use, and a package-managed pool is stopped before the next phase starts. Stage 2 fragment aggregation uses the same cap; do not set `metacell_args$fragment_args$workers` separately.

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
    min_cells = 500L,
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L
    )
  ),
  workers = workers
)
```

`pando_args$min_cells` defaults to `500L`, but it is not fixed. Any positive integer supplied by the user is retained and becomes the Stage 1 filtering threshold.

A cell type with at least two retained conditions uses the common-dictionary condition GRN. RegCompass initializes a separate Pando object for each broad cell type, parallelizes the pooled-background and per-condition candidate-discovery jobs, waits at a strict barrier, unions exact `(target, TF, region)` triples into one frozen dictionary for that cell type, and only then parallelizes the condition × cell-type fixed-dictionary Gaussian identity GLMs. Different cell types never share or merge Pando peak/motif feature spaces. The main workflow uses `tf_cor = 0.1` and `peak_cor = 0.05` for candidate discovery. Final condition edges are active when their target fit has `fit_status == "ok"`, the coefficient is estimable, and BH-adjusted `padj < 0.05`; no second post-fit coefficient-size, correlation, or model-R² gate is applied.

The fitted Pando `estimate` remains in Stage 1 outputs for statistical inference, direction and audit, but its absolute magnitude is no longer a Layer-1 penalty weight. For every active edge, Layer 1 uses only `sign(estimate)`. The paired-cell predictor `TF * ATAC` is robustly self-scaled across all fitted cells of the same broad cell type and bounded with `tanh`. Signed edge activities for each target/condition are averaged across that target's active Pando edges before the exact Stage-2 SuperCell membership aggregation is finalized. This prevents target degree or a condition-specific excess of significant edges from increasing the regulatory score by itself, while retaining paired-cell TF/ATAC co-occurrence and condition-pure membership. No new public parameter is introduced.

A cell type with one retained condition uses standard Pando. Standard Pando is parallelized across broad cell types, and each individual Pando fit is kept single-process to avoid nested oversubscription. The `tf_cor` and `peak_cor` values are passed to `Pando::infer_grn()` exactly as routed from `pando_infer_args`; the defaults are fixed at `tf_cor = 0.1` and `peak_cor = 0.05`, with no cell-count-dependent correlation threshold. The standard-Pando Layer-1 route uses the same sign-only, self-scaled and active-edge-degree-normalized paired-cell projection, so coefficient magnitude is not reintroduced when a cell type has only one effective condition.

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
    min_cells_per_stratum = 500L,
    min_metacell_size = 10L,
    min_metacells_per_stratum = 2L
  ),
  workers = workers
)
```

`metacell_args$min_cells_per_stratum` also defaults to `500L` and is user-configurable. If Stage 1 is intentionally run with a lower `min_cells`, set the Stage 2 threshold explicitly as needed rather than relying on the `500L` default.

One WNN graph is constructed per broad cell type. Final metacells remain condition-pure. If raw fragment files are supplied, `SuperCell::AggregateFragmentFile()` receives the same protected top-level worker cap.

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
  gene_half_saturation = 1,
  workers = workers
)
```

After the degree-controlled signed paired-cell activities have been aggregated to metacells, the existing target-level robust projection calibration and bounded regulatory modifier are unchanged. Target RNA remains the baseline gene support, regulatory evidence remains a bounded odds correction, and GPR aggregation, COMPASS reaction costs and Layer-2 LP mathematics are unchanged.

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

With `model_mode = "meta_module_gem"`, omitting `model_completion` selects original MATLAB CORDA2 semantics.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  workers = workers,
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

Do not set `completion_time_limit` for CORDA2. The canonical CORDA2 option object and Layer 2 completion context both use `Inf`; a supplied structural time limit is rejected rather than silently changing original CORDA2 execution.

### CORDA2 step-level parallelism

The mathematical order remains fixed:

```text
Step 1 HC dependencies
  -> barrier / deterministic ordered reduce / state update
  -> release worker pool and target-local HiGHS engines / GC
Step 2.1 MC-NC dependencies
  -> barrier / deterministic ordered reduce / state update
  -> release pool / GC
Step 2.2 MC feasibility
  -> barrier / deterministic ordered reduce / state update
  -> release pool / GC
Step 3 HC-OT dependencies
  -> barrier / deterministic ordered reduce / state update
  -> release pool / GC
```

Within the current step, directional targets are parallelized up to the protected worker cap. Each directional target starts from a fresh solver engine so its incoming simplex basis cannot depend on the worker/chunk assignment; repeated solves belonging to that same target still reuse its engine. HiGHS remains one thread per worker. Results are reduced in the original directional order only after all targets in the step finish.

The console and Layer 2 progress files report the current CORDA2 step, completed targets, total targets and `remaining=` count. Worker pools are short-lived by step; a new pool is created only after the preceding step has fully completed.

## 7. Result assembly

Use the Layer 2 output directly or continue with the one-shot result assembly functions. Targeted post-analysis can reuse either audited CORDA2 or supplementary FASTCORE Stage 5 union GEMs without structural reconstruction; `rc_regcompass_step_target_union()` uses the same top-level `workers` contract.
