# Tutorial 2: restartable workflow

Each stage writes a checkpoint to its output directory. Reuse the same input object, GEM, metadata columns, assays, and medium definition when restarting downstream stages. The code below shows the user-adjustable parameters used by the current workflow; fixed internal contracts are not exposed here.

```r
workers <- 10L
```

`workers` is the single workflow-level process cap. RegCompass selects the platform-specific backend automatically. For condition GRNs, the cap is divided across concurrently running cell types and each Pando fit can use its assigned remainder for target-level work; completed target pools are released immediately.

## 1. Regulatory evidence

```r
step1 <- rc_regcompass_step_grn(
  object = A,
  gem = gem,
  outdir = "run/01_grn",
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = NULL,
  species = "auto",
  condition_col = "condition",
  celltype_col = "cell_type",
  cell_type = NULL,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  pando_args = list(
    min_cells = 500L,
    pando_initiate_args = list(
      exclude_exons = TRUE
    ),
    pando_motif_args = list(
      reuse_cache = TRUE
    ),
    pando_infer_args = list(
      tf_cor = 0.1,
      peak_cor = 0.05,
      adjust_method = "BH",
      padj_threshold = 0.05,
      rank_action = "mark",
      min_residual_df = 1L,
      condition_ridge_control = list(
        lambda_grid = 10^seq(-3, 2, length.out = 9L),
        lambda_rule = "1se",
        fusion_ratio = 1,
        cv_folds = 5L,
        seed = 1L,
        scale_floor = 1e-8
      ),
      method = "glm"
    )
  ),
  workers = workers,
  progress = TRUE
)
```

For broad cell types with at least two retained conditions, RegCompass always uses the condition-comparable multi-task ridge path; `condition_ridge_control` controls its CV and shrinkage. `method` applies only to the standard-Pando route used by cell types with one effective condition. The default standard method remains `"glm"`.

To use the same ridge numerical backend for standard Pando, set `method = "ridge"` and optionally add `ridge_control`:

```r
pando_args_ridge <- list(
  min_cells = 500L,
  pando_infer_args = list(
    tf_cor = 0.1,
    peak_cor = 0.05,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = "mark",
    min_residual_df = 1L,
    method = "ridge",
    ridge_control = list(
      lambda_grid = 10^seq(-3, 2, length.out = 9L),
      lambda_rule = "1se",
      fusion_ratio = 1,
      cv_folds = 5L,
      seed = 1L,
      scale_floor = 1e-8
    ),
    condition_ridge_control = list(
      lambda_grid = 10^seq(-3, 2, length.out = 9L),
      lambda_rule = "1se",
      fusion_ratio = 1,
      cv_folds = 5L,
      seed = 1L,
      scale_floor = 1e-8
    )
  )
)
```

For standard ridge, `fusion_ratio` is accepted for a shared control schema but has no numerical effect because the standard model has one task. Standard-only Pando controls such as `peak_to_gene_method`, `upstream`, `downstream`, `extend`, `only_tss`, `alpha`, `family`, and other backend-specific arguments may also be supplied through `pando_infer_args` when that standard backend uses them.

## 2. Multimodal metacells

```r
step2 <- rc_regcompass_step_metacells(
  object = A,
  grn = step1,
  outdir = "run/02_metacells",
  condition_col = "condition",
  celltype_col = "cell_type",
  cell_type = NULL,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  fragment_files = NULL,
  metacell_args = list(
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 500L,
    min_metacell_size = 1L,
    min_metacells_per_stratum = 1L,
    k.knn = 30L,
    kith = NULL,
    kernel = TRUE,
    graph.name = NULL,
    metacellNormalization = FALSE,
    avg.in.data = FALSE,
    verbose = FALSE,
    fragment_args = list(
      rows_per_chunk = 10000000L,
      bgzip_path = NULL,
      tabix_path = NULL,
      process_n = 2000L,
      call_peaks = TRUE,
      macs2_path = NULL,
      effective_genome_size = NULL,
      peak_calling_args = list()
    )
  ),
  workers = workers,
  progress = TRUE
)
```

`fragment_args` is relevant only when `fragment_files` is supplied. Fragment worker count is deliberately not a separate parameter; the top-level `workers` value remains the single cap.

## 3. Reaction meta-modules

```r
step3 <- rc_regcompass_step_meta_modules(
  grn = step1,
  metacells = step2,
  gem = gem,
  outdir = "run/03_meta_modules",
  meta_module_args = list(
    subsystem_table = NULL
  ),
  progress = TRUE
)
```

`subsystem_table` is the current optional Stage 3 customization. The expansion order and GPR-completeness rules are fixed workflow contracts.

## 4. Layer 1 reaction evidence

```r
step4 <- rc_regcompass_step_layer1(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  gem = gem,
  outdir = "run/04_layer1",
  gpr_and_method = "min",
  gene_half_saturation = 1,
  workers = workers,
  progress = TRUE
)
```

`gpr_and_method` accepts `"min"`, `"median"`, or `"mean"`. `gene_half_saturation` is user-adjustable.

## 5. Medium scenarios

```r
medium_scenarios <- rc_make_medium_scenarios(
  gem = gem,
  scenario = "normal_human_plasma",
  species = "human",
  custom_medium = NULL,
  custom_metabolites = NULL,
  uptake_scale = 1,
  exchange_roles = "exchange",
  condition = "all",
  exchange_limit = 1,
  strict_preset_matching = TRUE
)
```

Built-in scenarios are `normal_human_plasma`, `mouse_plasma`, `high_glucose`, `low_glucose`, `high_lactate`, `low_lactate`, `low_glutamine`, and `custom`. See [medium-presets.md](medium-presets.md) for composition and provenance.

## 6. Layer 2 structural model and directional scoring

The default `meta_module_gem` route uses CORDA2. Its adjustable parameters are shown explicitly below; CORDA2 has no completion time limit.

```r
step5 <- rc_regcompass_step_layer2(
  layer1 = step4,
  meta_modules = step3,
  gem = gem,
  medium_scenarios = medium_scenarios,
  outdir = "run/05_layer2",
  model_mode = "meta_module_gem",
  layer2_args = list(
    omega = 0.95,
    target_direction = "both",
    solver = "highs",
    flux_threshold = 1e-8,
    model_params = list(
      model_completion = "corda2",
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
  ),
  workers = workers,
  progress = TRUE
)
```

`target_direction` accepts `"both"`, `"forward"`, or `"reverse"`; `solver` accepts the installed supported LP backends. The supplementary FASTCORE route is selected with `model_completion = "fastcore"`; its adjustable structural controls are `completion_time_limit`, `fastcore_epsilon`, `max_support_reactions`, and `strict`. See [layer2-corda.md](layer2-corda.md) for the CORDA2 parameter contract.

## 7. Result assembly

```r
result <- rc_regcompass_step_results(
  grn = step1,
  metacells = step2,
  meta_modules = step3,
  layer1 = step4,
  layer2 = step5,
  gem = gem,
  outdir = "run/06_results",
  species = "human",
  progress = TRUE
)
```

The stage outputs are restartable R objects. Post-analysis functions are listed in [tutorial-04-post-analysis.md](tutorial-04-post-analysis.md) and [functions.md](functions.md). Equations and quantitative definitions are maintained only in [mathematical-model.md](mathematical-model.md).
