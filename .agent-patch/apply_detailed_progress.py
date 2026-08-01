from pathlib import Path


def replace_once(path, old, new):
    file = Path(path)
    text = file.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one match, found {count}: {old[:100]!r}")
    file.write_text(text.replace(old, new, 1))


def insert_after_once(path, anchor, addition):
    replace_once(path, anchor, anchor + addition)


# Stage-1 orchestration: expose deterministic high-level phases and pass the
# monitor into the canonical condition-GRN bridge. The public API is unchanged.
replace_once(
    "R/stepwise_workflow.R",
    '  monitor <- .rc_step_monitor_start("grn", outdir, progress)\n',
    '  monitor <- .rc_step_monitor_start(\n'
    '    "grn", outdir, progress, total_parts = 12L\n'
    '  )\n'
)
insert_after_once(
    "R/stepwise_workflow.R",
    '  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)\n',
    '  .rc_step_monitor_event(\n'
    '    monitor, "input_validation",\n'
    '    "validating Stage 1 arguments and GEM", current = 1L\n'
    '  )\n'
)
insert_after_once(
    "R/stepwise_workflow.R",
    '  species <- .rc_infer_gem_species(gem, species)\n  rc_validate_gem(gem)\n',
    '  .rc_step_monitor_event(\n'
    '    monitor, "gem_validated", "GEM contract validated", current = 1L,\n'
    '    context = list(species = species)\n'
    '  )\n'
)
insert_after_once(
    "R/stepwise_workflow.R",
    '  object <- design$object\n  effective_condition_col <- design$condition_col\n',
    '  .rc_step_monitor_event(\n'
    '    monitor, "design_resolution",\n'
    '    "resolved condition-aware versus standard Pando route", current = 2L,\n'
    '    context = list(\n'
    '      analysis_mode = design$analysis_mode,\n'
    '      condition_col = effective_condition_col %||% "<none>",\n'
    '      condition_levels = paste(design$condition_levels, collapse = ",")\n'
    '    )\n'
    '  )\n'
)
insert_after_once(
    "R/stepwise_workflow.R",
    '  object <- .rc_normalize_single_cell_grn_object(\n'
    '    object,\n'
    '    condition_col = effective_condition_col,\n'
    '    celltype_col = celltype_col,\n'
    '    rna_assay = rna_assay,\n'
    '    atac_assay = atac_assay\n'
    '  )\n',
    '  .rc_step_monitor_event(\n'
    '    monitor, "single_cell_normalization",\n'
    '    "RNA and ATAC inputs normalized for GRN inference", current = 3L,\n'
    '    context = list(\n'
    '      cells = ncol(object),\n'
    '      cell_types = length(unique(as.character(\n'
    '        object@meta.data[[celltype_col]]\n'
    '      ))),\n'
    '      rna_assay = rna_assay,\n'
    '      atac_assay = atac_assay\n'
    '    )\n'
    '  )\n'
)
replace_once(
    "R/stepwise_workflow.R",
    '    infer_args$candidate_screen <- infer_args$candidate_screen %||% "motif_domain"\n'
    '    infer_args$parallel <- FALSE\n'
    '    call_args$pando_infer_args <- infer_args\n'
    '    call_args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE\n'
    '    grn_result <- do.call(.rc_fit_condition_grns_by_cell_type, call_args)\n',
    '    infer_args$candidate_screen <- infer_args$candidate_screen %||% "motif_domain"\n'
    '    infer_args$parallel <- FALSE\n'
    '    infer_args$verbose <- infer_args$verbose %||%\n'
    '      .rc_progress_enabled(progress)\n'
    '    call_args$pando_infer_args <- infer_args\n'
    '    call_args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE\n'
    '    call_args$progress_monitor <- monitor\n'
    '    .rc_step_monitor_event(\n'
    '      monitor, "pando_configuration",\n'
    '      "configured fused C++ condition-GRN runtime", current = 4L,\n'
    '      context = list(\n'
    '        candidate_screen = infer_args$candidate_screen,\n'
    '        outer_folds = infer_args$outer_nfolds %||% 5L,\n'
    '        inner_folds = infer_args$inner_nfolds %||% 5L,\n'
    '        nlambda = infer_args$nlambda %||% 50L,\n'
    '        lambda_selection = infer_args$lambda_selection %||% "lambda.1se",\n'
    '        backend = "cpp_hybrid_gram_sufficient_statistics"\n'
    '      )\n'
    '    )\n'
    '    grn_result <- do.call(.rc_fit_condition_grns_by_cell_type, call_args)\n'
)
replace_once(
    "R/stepwise_workflow.R",
    '  } else {\n'
    '    call_args$pando_infer_args <- infer_args\n'
    '    call_args$parallel <- isTRUE(parallel)\n'
    '    grn_result <- do.call(.rc_fit_standard_pando_by_cell_type, call_args)\n'
    '  }\n',
    '  } else {\n'
    '    infer_args$verbose <- infer_args$verbose %||%\n'
    '      .rc_progress_enabled(progress)\n'
    '    call_args$pando_infer_args <- infer_args\n'
    '    call_args$parallel <- isTRUE(parallel)\n'
    '    .rc_step_monitor_event(\n'
    '      monitor, "standard_pando",\n'
    '      "dispatching original Pando infer_grn workflow", current = 5L\n'
    '    )\n'
    '    grn_result <- do.call(.rc_fit_standard_pando_by_cell_type, call_args)\n'
    '    .rc_step_monitor_event(\n'
    '      monitor, "standard_pando_complete",\n'
    '      "original Pando workflow completed", current = 10L\n'
    '    )\n'
    '  }\n'
)
insert_after_once(
    "R/stepwise_workflow.R",
    '  grn_result$atac_assay <- atac_assay\n',
    '  .rc_step_monitor_event(\n'
    '    monitor, "stage_contract",\n'
    '    "assembled RegCompass Stage 1 GRN contract", current = 11L,\n'
    '    context = list(analysis_mode = design$analysis_mode)\n'
    '  )\n'
)

# Condition-GRN bridge: report each biologically meaningful preprocessing and
# numerical phase while preserving the original function/API semantics.
replace_once(
    "R/condition_grn_contract.R",
    '    save_pando_objects = TRUE, BPPARAM = NULL,\n'
    '    species = c("auto", "human", "mouse")) {\n',
    '    save_pando_objects = TRUE, BPPARAM = NULL,\n'
    '    progress_monitor = NULL,\n'
    '    species = c("auto", "human", "mouse")) {\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  species <- .rc_infer_gem_species(gem, species)\n  rc_validate_gem(gem)\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "pando_contract_check",\n'
    '    "validating Pando runtime and paired cell metadata", current = 5L,\n'
    '    context = list(species = species)\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  if (!identical(pando_infer_args$candidate_screen, "motif_domain") ||\n'
    '      !identical(pando_infer_args$condition_weight, "equal") ||\n'
    '      !isTRUE(pando_infer_args$scale)) {\n'
    '    stop(\n'
    '      "Condition comparability requires candidate_screen=\'motif_domain\', condition_weight=\'equal\', and scale=TRUE.",\n'
    '      call. = FALSE\n'
    '    )\n'
    '  }\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "condition_design",\n'
    '    "condition-comparable fit controls validated", current = 5L,\n'
    '    context = list(\n'
    '      conditions = length(unique(as.character(\n'
    '        object@meta.data[[condition_col]]\n'
    '      ))),\n'
    '      cell_types = length(unique(as.character(\n'
    '        object@meta.data[[celltype_col]]\n'
    '      ))),\n'
    '      outer_folds = pando_infer_args$outer_nfolds,\n'
    '      inner_folds = pando_infer_args$inner_nfolds,\n'
    '      nlambda = pando_infer_args$nlambda,\n'
    '      lambda_selection = pando_infer_args$lambda_selection\n'
    '    )\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  if (!length(target_genes)) {\n'
    '    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)\n'
    '  }\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "target_selection",\n'
    '    "resolved metabolic target genes shared by RNA and GEM", current = 6L,\n'
    '    context = list(targets = length(target_genes))\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  object <- filtered$object\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "atac_feature_filter",\n'
    '    "removed globally zero ATAC features", current = 6L,\n'
    '    context = list(\n'
    '      cells = ncol(object),\n'
    '      removed_atac_features = filtered$n_removed %||% NA_integer_\n'
    '    )\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "candidate_initialization",\n'
    '    "initialized Pando regulatory candidate space", current = 7L,\n'
    '    context = list(targets = length(target_genes))\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  grn <- do.call(Pando::find_motifs, c(motif, pando_motif_args))\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "motif_mapping",\n'
    '    "completed motif-to-peak and TF mapping", current = 8L\n'
    '  )\n'
)
replace_once(
    "R/condition_grn_contract.R",
    '  infer[names(pando_infer_args)] <- NULL\n'
    '  grn <- do.call(Pando::infer_condition_grn, c(infer, pando_infer_args))\n'
    '  extracted <- .rc_extract_condition_grn_contract(\n',
    '  infer[names(pando_infer_args)] <- NULL\n'
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "nested_cv",\n'
    '    "running fused per-target outer/inner CV, path, refit and validation",\n'
    '    current = 9L,\n'
    '    context = list(\n'
    '      targets = length(target_genes),\n'
    '      cell_types = length(unique(as.character(\n'
    '        object@meta.data[[celltype_col]]\n'
    '      ))),\n'
    '      conditions = length(unique(as.character(\n'
    '        object@meta.data[[condition_col]]\n'
    '      ))),\n'
    '      outer_folds = pando_infer_args$outer_nfolds,\n'
    '      inner_folds = pando_infer_args$inner_nfolds,\n'
    '      nlambda = pando_infer_args$nlambda,\n'
    '      solver = "hybrid_gram_or_sparse_matrix_free",\n'
    '      validation = "exact_sufficient_statistics",\n'
    '      oof = "outer_selected_model_only"\n'
    '    )\n'
    '  )\n'
    '  grn <- do.call(Pando::infer_condition_grn, c(infer, pando_infer_args))\n'
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "nested_cv_complete",\n'
    '    "Pando fused target engine completed", current = 10L\n'
    '  )\n'
    '  extracted <- .rc_extract_condition_grn_contract(\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '    min_model_rsq = min_model_rsq\n  )\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "contract_extraction",\n'
    '    "validated and extracted ConditionGRNFit contracts", current = 10L,\n'
    '    context = list(\n'
    '      fitted_cell_types = length(extracted$fit_contracts),\n'
    '      all_edges = nrow(extracted$condition_all),\n'
    '      active_edges = nrow(extracted$condition_active)\n'
    '    )\n'
    '  )\n'
)
insert_after_once(
    "R/condition_grn_contract.R",
    '  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)\n',
    '  .rc_step_monitor_event(\n'
    '    progress_monitor, "grn_artifacts",\n'
    '    "writing Stage 1 GRN tables and fit contracts", current = 11L,\n'
    '    context = list(groups = nrow(status))\n'
    '  )\n'
)

# Documentation and regression coverage.
insert_after_once(
    "NEWS.md",
    "# RegCompassR 2.2.1\n",
    "\n- Added structured stage progress for Stage 1. Console events now report "
    "phase, percent, elapsed time, cell types, conditions, metabolic targets, "
    "outer/inner folds, lambda path size, solver route, refit/validation and "
    "OOF status. Every run writes `step_progress.tsv`, including when console "
    "messages are disabled.\n"
    "- Stage 1 now enables Pando target-level verbose output when RegCompass "
    "progress is enabled, while retaining the existing public API and "
    "fail-fast behavior.\n"
)
insert_after_once(
    "docs/tutorial-02-stepwise-audit.md",
    "[functions.md](functions.md).\n",
    "\n## Detailed progress and audit log\n\n"
    "All public stages accept `progress = TRUE`. Stage 1 reports input "
    "validation, design resolution, normalization, Pando runtime checks, "
    "target selection, ATAC filtering, candidate initialization, motif "
    "mapping, nested CV, contract extraction and artifact writing. The "
    "long-running nested-CV event includes the number of cell types, "
    "conditions, metabolic targets, outer/inner folds and lambda values.\n\n"
    "```r\n"
    "options(RegCompassR.progress = TRUE)\n"
    "step1 <- rc_regcompass_step_grn(..., progress = TRUE)\n"
    "progress_log <- read.delim(\n"
    "  \"RegCompass_steps/01_grn/step_progress.tsv\",\n"
    "  check.names = FALSE\n"
    ")\n"
    "progress_log[, c(\"phase\", \"elapsed_hms\", \"detail\", \"context\")]\n"
    "```\n\n"
    "`step_progress.tsv` is written even when `progress = FALSE`; that setting "
    "only suppresses console messages. Pando target-level messages are enabled "
    "automatically when Stage 1 progress is enabled. Errors terminate "
    "immediately and the last audit row is `stage_error`.\n"
)
insert_after_once(
    "vignettes/regcompass-workflow.Rmd",
    "support and directional metabolic scores.\n",
    "\nEach stage writes `step_timing.tsv` and a structured `step_progress.tsv`. "
    "With `progress = TRUE`, the console also reports phase, percent and "
    "elapsed time. Stage 1 additionally records cell types, conditions, "
    "metabolic targets, outer/inner folds, lambda path size, hybrid solver "
    "route, contract extraction and artifact writing. Setting `progress = "
    "FALSE` suppresses console output but retains the audit log.\n"
)

Path("tests/testthat/test-detailed-stage-progress.R").write_text(r'''test_that("stage progress reports percent, elapsed time and structured context", {
  outdir <- tempfile("regcompass-progress-")
  monitor <- .rc_step_monitor_start(
    "grn", outdir, progress = TRUE, total_parts = 4L
  )
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)

  messages <- testthat::capture_messages(
    .rc_step_monitor_event(
      monitor,
      phase = "nested_cv",
      detail = "running fused target engine",
      current = 2L,
      context = list(
        cell_types = 3L,
        conditions = 2L,
        targets = 120L,
        outer_folds = 5L,
        inner_folds = 5L,
        nlambda = 50L
      )
    )
  )
  text <- paste(messages, collapse = "\n")
  expect_match(text, "phase=nested_cv", fixed = TRUE)
  expect_match(text, "50.0%", fixed = TRUE)
  expect_match(text, "elapsed=", fixed = TRUE)
  expect_match(text, "targets=120", fixed = TRUE)

  log <- utils::read.delim(
    file.path(outdir, "step_progress.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_true(all(c(
    "sequence", "timestamp", "stage", "phase", "status", "current",
    "total", "percent", "elapsed_seconds", "elapsed_hms", "detail",
    "context"
  ) %in% colnames(log)))
  expect_true(all(c("stage_start", "nested_cv") %in% log$phase))
  expect_match(log$context[log$phase == "nested_cv"], "outer_folds=5")
})

test_that("progress FALSE suppresses messages but retains the audit log", {
  outdir <- tempfile("regcompass-progress-silent-")
  expect_silent({
    monitor <- .rc_step_monitor_start(
      "grn", outdir, progress = FALSE, total_parts = 3L
    )
    .rc_step_monitor_event(
      monitor, "motif_mapping", "mapped motifs", current = 2L,
      context = list(targets = 10L)
    )
    .rc_step_monitor_fail(monitor)
  })
  log <- utils::read.delim(
    file.path(outdir, "step_progress.tsv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_true(all(c("stage_start", "motif_mapping", "stage_error") %in%
    log$phase))
})

test_that("Stage 1 source forwards monitor and verbose status to Pando", {
  source_text <- paste(
    readLines(test_path("..", "..", "R", "stepwise_workflow.R"), warn = FALSE),
    collapse = "\n"
  )
  bridge_text <- paste(
    readLines(test_path("..", "..", "R", "condition_grn_contract.R"),
      warn = FALSE),
    collapse = "\n"
  )
  expect_match(source_text, "call_args\\$progress_monitor <- monitor")
  expect_match(source_text, "infer_args\\$verbose <- infer_args\\$verbose")
  expect_match(bridge_text, "phase = \\"nested_cv\\"")
  expect_match(bridge_text, "exact_sufficient_statistics")
  expect_match(bridge_text, "outer_selected_model_only")
})
''')
