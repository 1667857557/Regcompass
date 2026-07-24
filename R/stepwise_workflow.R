.rc_workflow_signature <- function(x) {
  params <- x$params %||% x$workflow_params %||% list()
  params[c("condition_col", "celltype_col", "rna_assay", "atac_assay")]
}

.rc_validate_grn_metacell_group_coverage <- function(
    grn_result, metacell_meta,
    condition_col = "condition", celltype_col = "cell_type") {
  group_cols <- c(condition_col, celltype_col)
  status <- grn_result$sample_status
  if (!is.data.frame(status) ||
      !all(c(group_cols, "status") %in% colnames(status))) {
    stop("GRN status is incomplete for condition-by-cell-type coverage validation.",
         call. = FALSE)
  }
  if (!is.data.frame(metacell_meta) ||
      !all(group_cols %in% colnames(metacell_meta))) {
    stop("Metacell metadata are incomplete for GRN coverage validation.",
         call. = FALSE)
  }
  status$.group_id <- rc_make_stratum_id(status, group_cols)
  metacell_meta$.group_id <- rc_make_stratum_id(metacell_meta, group_cols)

  grn_rows <- split(seq_len(nrow(status)), status$.group_id)
  grn_summary <- do.call(rbind, lapply(grn_rows, function(rows) {
    one <- status[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      grn_status = paste(sort(unique(as.character(one$status))), collapse = ";"),
      n_single_cells = sum(as.numeric(one$n_cells %||% 0), na.rm = TRUE),
      n_significant_edges = sum(
        as.numeric(one$n_significant_edges %||% 0), na.rm = TRUE
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  metacell_rows <- split(
    seq_len(nrow(metacell_meta)), metacell_meta$.group_id
  )
  metacell_summary <- do.call(rbind, lapply(metacell_rows, function(rows) {
    one <- metacell_meta[rows, , drop = FALSE]
    values <- one[1L, group_cols, drop = FALSE]
    purity <- if ("dominant_celltype_fraction" %in% colnames(one)) {
      suppressWarnings(as.numeric(one$dominant_celltype_fraction))
    } else {
      rep(NA_real_, nrow(one))
    }
    mixed <- if ("mixed_celltype_metacell" %in% colnames(one)) {
      one$mixed_celltype_metacell %in% TRUE
    } else {
      rep(FALSE, nrow(one))
    }
    data.frame(
      values,
      group_id = as.character(one$.group_id[[1L]]),
      n_metacells = nrow(one),
      median_dominant_celltype_fraction = if (all(is.na(purity))) {
        NA_real_
      } else {
        stats::median(purity, na.rm = TRUE)
      },
      min_dominant_celltype_fraction = if (all(is.na(purity))) {
        NA_real_
      } else {
        min(purity, na.rm = TRUE)
      },
      n_mixed_celltype_metacells = sum(mixed, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))

  coverage <- merge(
    grn_summary, metacell_summary,
    by = c(group_cols, "group_id"), all = TRUE, sort = TRUE
  )
  coverage$grn_available <- !is.na(coverage$grn_status) &
    coverage$grn_status == "ok" & coverage$n_significant_edges > 0
  coverage$metacells_available <- !is.na(coverage$n_metacells) &
    coverage$n_metacells > 0
  coverage$coverage_complete <- coverage$grn_available &
    coverage$metacells_available
  invalid <- coverage[!coverage$coverage_complete, , drop = FALSE]
  if (nrow(invalid)) {
    stop(
      "GRN and metacell condition-by-cell-type groups do not align: ",
      paste(invalid$group_id, collapse = "; "),
      ". Every scored metacell group requires a successful GRN with significant edges, and every GRN group requires at least one metacell.",
      call. = FALSE
    )
  }
  rownames(coverage) <- NULL
  coverage
}

#' Infer condition-by-cell-type Pando GRNs from single cells
#' @export
rc_regcompass_step_grn <- function(
    object, gem, outdir, pfm, genome,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL) {
  if (!is.list(pando_args)) stop("`pando_args` must be a list.", call. = FALSE)
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  object <- .rc_normalize_single_cell_grn_object(
    object, condition_col = condition_col, celltype_col = celltype_col,
    rna_assay = rna_assay, atac_assay = atac_assay
  )
  reserved <- intersect(names(pando_args), c(
    "object", "gem", "outdir", "pfm", "genome", "condition_col",
    "celltype_col", "rna_assay", "atac_assay", "BPPARAM"
  ))
  if (length(reserved)) {
    stop("`pando_args` cannot override workflow fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  defaults <- list(
    object = object, gem = gem, outdir = outdir,
    pfm = pfm, genome = genome, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay,
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE,
    on_group_error = "stop"
  )
  defaults[names(pando_args)] <- NULL
  grn_result <- do.call(
    .rc_run_condition_single_cell_grns, c(defaults, pando_args)
  )
  answer <- list(
    grn_result = grn_result,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      condition_col = condition_col, celltype_col = celltype_col,
      rna_assay = rna_assay, atac_assay = atac_assay,
      pando_args = pando_args, parallel = parallel
    )
  )
  class(answer) <- c("regcompass_grn_step", "list")
  saveRDS(answer, file.path(outdir, "step_grn.rds"))
  answer
}

#' Build condition-only SuperCell2 metacells
#' @export
rc_regcompass_step_metacells <- function(
    object, outdir,
    sample_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (identical(fragment_files, FALSE) || is.null(fragment_files)) {
    object <- .rc_clear_signac_fragments(object, atac_assay = atac_assay)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pooled <- .rc_make_condition_pooled_metacells(
    object = object, outdir = outdir, sample_col = sample_col,
    condition_col = condition_col, celltype_col = celltype_col,
    rna_assay = rna_assay, atac_assay = atac_assay,
    fragment_files = fragment_files, metacell_args = metacell_args
  )
  metacell_object <- .rc_normalize_condition_metacell_object(
    pooled, rna_assay, atac_assay
  )
  if (!setequal(
    colnames(metacell_object), as.character(pooled$metacell_meta$metacell_id)
  )) {
    stop("Merged metacell object and metadata contain different units.",
         call. = FALSE)
  }
  .rc_write_tsv_gz(
    pooled$metacell_meta,
    file.path(outdir, "metacell_metadata.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$membership,
    file.path(outdir, "metacell_membership.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition,
    file.path(outdir, "metacell_celltype_composition.tsv.gz")
  )
  .rc_write_tsv_gz(
    pooled$celltype_composition_summary,
    file.path(outdir, "metacell_celltype_summary.tsv.gz")
  )
  saveRDS(metacell_object, file.path(outdir, "merged_metacell_object.rds"))
  answer <- list(
    pooled = pooled,
    metacell_object = metacell_object,
    params = list(
      input_sample_col = sample_col,
      sample_col = pooled$analysis_sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = modifyList(list(gamma = 75L), metacell_args)
    )
  )
  class(answer) <- c("regcompass_metacell_step", "list")
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}

#' Construct core reactions and condition-specific meta-modules from GRNs
#' @export
rc_regcompass_step_meta_modules <- function(
    grn, metacells, gem, outdir,
    layer1_args = list()) {
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  if (!identical(.rc_workflow_signature(grn),
                 .rc_workflow_signature(metacells))) {
    stop("GRN and metacell stages use different metadata or assay settings.",
         call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  validated_gem <- rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_coverage <- .rc_validate_grn_metacell_group_coverage(
    grn_result = grn$grn_result,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = metacells$params$condition_col,
    celltype_col = metacells$params$celltype_col
  )
  .rc_write_tsv_gz(
    group_coverage,
    file.path(outdir, "grn_metacell_group_coverage.tsv.gz")
  )
  condition_modules <- .rc_build_condition_meta_modules(
    grn$grn_result, gem, outdir, layer1_args
  )
  condition_modules$grn_metacell_group_coverage <- group_coverage
  if (!is.data.frame(condition_modules$reaction_membership) ||
      !nrow(condition_modules$reaction_membership)) {
    stop("Meta-module construction produced no reaction membership.",
         call. = FALSE)
  }
  missing <- setdiff(
    unique(as.character(condition_modules$reaction_membership$reaction_id)),
    colnames(validated_gem$S)
  )
  if (length(missing)) {
    stop("Meta-module reactions absent from the GEM: ",
         paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
  }
  global_modules <- .rc_merge_stratum_meta_modules(list(list(
    group_id = "condition_pooled", grn_meta_modules = condition_modules
  )))
  if (!is.data.frame(global_modules$global_core_reactions) ||
      !nrow(global_modules$global_core_reactions)) {
    stop("No complete-GPR global core reactions remain after module merging.",
         call. = FALSE)
  }
  answer <- list(
    condition_modules = condition_modules,
    global_modules = global_modules,
    group_coverage = group_coverage,
    workflow_params = metacells$params,
    grn_params = grn$params,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(layer1_args = layer1_args)
  )
  class(answer) <- c("regcompass_meta_module_step", "list")
  saveRDS(condition_modules, file.path(outdir, "condition_meta_modules.rds"))
  saveRDS(global_modules, file.path(outdir, "global_meta_modules.rds"))
  saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  answer
}

#' Build integrated RNA+ATAC reaction expression
#' @export
rc_regcompass_step_layer1 <- function(
    metacells, meta_modules, gem, outdir,
    regulatory_alpha = 1,
    tau = 0.20,
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE,
    BPPARAM = NULL) {
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop("Metacell and meta-module stages use different workflow settings.",
         call. = FALSE)
  }
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  params <- metacells$params
  layer1 <- .rc_build_condition_pooled_layer1(
    metacell_object = metacells$metacell_object,
    meta_modules = meta_modules$condition_modules,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    sample_col = params$sample_col,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    atac_assay = params$atac_assay,
    regulatory_alpha = regulatory_alpha,
    gpr_tau = tau,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(layer1) <- c("regcompass_layer1_step", "list")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}

#' Run directional COMPASS-like LP scoring
#' @export
rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {
  model_mode <- match.arg(model_mode)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!is.list(layer2_args)) stop("`layer2_args` must be a list.", call. = FALSE)
  params <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  layer2_args$model_params$cache_dir <- file.path(
    outdir, "model_cache", model_mode
  )
  reserved <- intersect(names(layer2_args), c(
    "layer1", "gem", "mode", "unit", "reaction_membership",
    "core_reactions", "target_reactions", "medium_scenarios",
    "sample_col", "condition_col", "celltype_col", "BPPARAM",
    "parallel", "penalty_weights"
  ))
  if (length(reserved)) {
    stop("`layer2_args` cannot override workflow fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  .rc_require_lp_solver(solver)
  global <- meta_modules$global_modules
  targets <- unique(as.character(global$global_core_reactions$reaction_id))
  missing_expression <- setdiff(targets, rownames(layer1$reaction_expression))
  if (length(missing_expression)) {
    stop("Global core reactions are absent from Layer 1 expression: ",
         paste(utils::head(missing_expression, 10L), collapse = ", "),
         call. = FALSE)
  }
  defaults <- list(
    layer1 = layer1, gem = gem,
    target_reactions = targets,
    medium_scenarios = medium_scenarios, mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      global$global_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      global$global_core_reactions
    } else {
      NULL
    },
    unit = "metacell", sample_col = params$sample_col,
    condition_col = params$condition_col, celltype_col = params$celltype_col,
    parallel = parallel, BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  answer <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(w) {
      if (grepl(
        "Metacell-level scores are descriptive pseudo-observations",
        conditionMessage(w), fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  answer$workflow_params <- params
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  answer$source_core_reactions <- global$global_core_reactions
  class(answer) <- c("regcompass_layer2_step", "list")
  .rc_validate_layer2_stage(
    answer, layer1 = layer1, workflow_params = params, gem = gem,
    argument = "layer2"
  )
  rc_export_microcompass(answer, outdir)
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}

#' Assemble final RegCompass results
#' @export
rc_regcompass_step_results <- function(
    grn, metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse")) {
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  params <- metacells$params
  if (!identical(params, meta_modules$workflow_params) ||
      !identical(.rc_workflow_signature(grn), .rc_workflow_signature(metacells))) {
    stop("Upstream stages use different workflow parameters.", call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = params, gem = gem,
    argument = "layer2"
  )
  species <- .rc_infer_gem_species(gem, species)
  comparison <- .rc_condition_penalty_comparison(
    layer2,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  conditions <- unique(as.character(
    metacells$pooled$metacell_meta[[params$condition_col]]
  ))
  completion <- .rc_feasibility_completion_metadata(layer2$model_mode)
  condition_fields <- intersect(c(
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "biological_reaction_membership", "reaction_membership",
    "local_completed_reaction_membership", "meta_module_summary",
    "local_fastcore_summary", "local_fastcore_diagnostics",
    "local_fastcore_completion_iterations", "local_fastcore_parent_scope",
    "analysis_group_unit", "grn_metacell_group_coverage"
  ), names(meta_modules$condition_modules))
  condition_modules <- meta_modules$condition_modules[condition_fields]
  result <- list(
    schema_version = "regcompass_grn_first_v2",
    version = "1.8.2", species = species, model_mode = layer2$model_mode,
    analysis_mode = comparison$analysis_mode,
    grn = grn$grn_result,
    metacells = metacells$pooled,
    layer1 = layer1,
    condition_grn_meta_modules = condition_modules,
    global_grn_meta_modules = meta_modules$global_modules,
    grn_meta_modules = meta_modules$global_modules,
    grn_metacell_group_coverage = meta_modules$group_coverage,
    microcompass = layer2,
    reaction_ranking = comparison$ranking,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    inference_policy = comparison$inference_policy,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      n_conditions = length(conditions),
      workflow_order = c(
        "single_cell_grn", "condition_metacells", "meta_modules",
        "layer1", "layer2"
      ),
      pando_grouping = c(params$condition_col, params$celltype_col),
      pando_peak_cor = grn$grn_result$normalization_policy$pando_peak_cor,
      metacell_grouping = params$condition_col,
      metacell_celltype_assignment =
        "supercell_label_guided_then_dominant_membership_audit",
      metacell_gamma = params$metacell_args$gamma,
      sample_weighting = "none",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      feasibility_completion = completion$feasibility_completion,
      feasibility_completion_stages = completion$feasibility_completion_stages,
      pando_normalization_policy = grn$grn_result$normalization_policy,
      penalty_formula = "1/(1+log2(1+E_multiome))",
      execution_mode = "stepwise"
    )
  )
  result <- .rc_ra_attach_to_result(
    result = result,
    gem = gem,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    result$reaction_catalog,
    file.path(outdir, "reaction_catalog.tsv.gz")
  )
  .rc_write_tsv_gz(
    result$reaction_evidence,
    file.path(outdir, "reaction_evidence_by_condition_celltype.tsv.gz")
  )
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
