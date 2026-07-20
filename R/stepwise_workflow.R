#' Build condition-pooled metacells as a standalone workflow step
#'
#' This is the first public stage of the inspectable RegCompass workflow. It
#' builds condition-by-cell-type SuperCell2 metacells, retains biological-sample
#' composition diagnostics, creates the normalized merged RNA+ATAC metacell
#' object, and persists the stage outputs before returning them.
#'
#' @param object A Seurat RNA+ATAC object.
#' @param outdir Output directory for this stage.
#' @param sample_col,condition_col,celltype_col Metadata columns.
#' @param rna_assay,atac_assay Assay names.
#' @param fragment_files Must be `FALSE` for the canonical condition-pooled path.
#' @param metacell_args Arguments passed to the SuperCell2 metacell builder.
#' @return A list with `pooled` diagnostics and `metacell_object`.
#' @export
rc_regcompass_step_metacells <- function(
    object, outdir,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  pooled <- withCallingHandlers(
    .rc_make_condition_pooled_metacells(
      object = object,
      outdir = outdir,
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = metacell_args,
      strict_biological_defaults = FALSE
    ),
    warning = function(warning) {
      if (grepl(
        "Condition-pooled analysis requires at least two biological samples",
        conditionMessage(warning),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
  metacell_object <- .rc_normalize_condition_metacell_object(
    pooled,
    rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  if (!setequal(
    colnames(metacell_object),
    as.character(pooled$metacell_meta$metacell_id)
  )) {
    stop(
      "Merged metacell object and pooled metadata contain different units.",
      call. = FALSE
    )
  }
  saveRDS(metacell_object, file.path(outdir, "merged_metacell_object.rds"))
  answer <- list(
    pooled = pooled,
    metacell_object = metacell_object,
    params = list(
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      fragment_files = fragment_files,
      metacell_args = metacell_args
    )
  )
  class(answer) <- c("regcompass_metacell_step", "list")
  saveRDS(answer, file.path(outdir, "step_metacells.rds"))
  answer
}

#' Build Pando GRNs and condition-specific meta-modules as a standalone step
#'
#' @param metacells Output from [rc_regcompass_step_metacells()].
#' @param gem A prepared GEM.
#' @param outdir Output directory for this stage.
#' @param pfm Motif position-frequency matrices.
#' @param genome Genome object matching the ATAC coordinates.
#' @param pando_args Arguments passed to Pando inference.
#' @param layer1_args Layer 1 structural arguments used by local FASTCORE.
#' @param BPPARAM BiocParallel parameter object or `FALSE`.
#' @return A list with condition-specific and merged global meta-modules.
#' @export
rc_regcompass_step_meta_modules <- function(
    metacells, gem, outdir, pfm, genome,
    pando_args = list(),
    layer1_args = list(),
    BPPARAM = FALSE) {
  if (!inherits(metacells, "regcompass_metacell_step")) {
    stop(
      "`metacells` must be the output of `rc_regcompass_step_metacells()`.",
      call. = FALSE
    )
  }
  params <- metacells$params
  validated_gem <- rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  condition_modules <- .rc_run_condition_pando_modules(
    metacell_object = metacells$metacell_object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    atac_assay = params$atac_assay,
    pando_args = pando_args,
    layer1_args = layer1_args,
    BPPARAM = BPPARAM
  )
  if (!is.data.frame(condition_modules$reaction_membership) ||
      !nrow(condition_modules$reaction_membership)) {
    stop(
      "Pando meta-module construction produced no reaction membership.",
      call. = FALSE
    )
  }
  missing_reactions <- setdiff(
    unique(as.character(condition_modules$reaction_membership$reaction_id)),
    colnames(validated_gem$S)
  )
  if (length(missing_reactions)) {
    stop(
      "Meta-module membership contains reactions absent from the GEM: ",
      paste(utils::head(missing_reactions, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  global_modules <- .rc_merge_stratum_meta_modules(list(list(
    group_id = "condition_pooled",
    grn_meta_modules = condition_modules
  )))
  if (!is.data.frame(global_modules$global_core_reactions) ||
      !nrow(global_modules$global_core_reactions)) {
    stop(
      "No complete-GPR global core reactions remain after module merging.",
      call. = FALSE
    )
  }
  answer <- list(
    condition_modules = condition_modules,
    global_modules = global_modules,
    params = list(
      pando_args = pando_args,
      layer1_args = layer1_args
    )
  )
  class(answer) <- c("regcompass_meta_module_step", "list")
  saveRDS(condition_modules, file.path(outdir, "condition_meta_modules.rds"))
  saveRDS(global_modules, file.path(outdir, "global_meta_modules.rds"))
  saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  answer
}

#' Build integrated RNA+ATAC reaction expression as a standalone step
#'
#' @param metacells Output from [rc_regcompass_step_metacells()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param gem A prepared GEM.
#' @param outdir Output directory for this stage.
#' @param regulatory_alpha Strength of the ATAC modifier on RNA support log-odds.
#' @param tau Boltzmann soft-min parameter for GPR AND rules.
#' @param gene_half_saturation Half-saturation for absolute RNA support.
#' @return The Layer 1 result containing RNA support, regulatory modifier,
#'   integrated gene support, and reaction expression.
#' @export
rc_regcompass_step_layer1 <- function(
    metacells, meta_modules, gem, outdir,
    regulatory_alpha = 1,
    tau = 0.20,
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1)) {
  if (!inherits(metacells, "regcompass_metacell_step")) {
    stop(
      "`metacells` must be the output of `rc_regcompass_step_metacells()`.",
      call. = FALSE
    )
  }
  if (!inherits(meta_modules, "regcompass_meta_module_step")) {
    stop(
      "`meta_modules` must be the output of `rc_regcompass_step_meta_modules()`.",
      call. = FALSE
    )
  }
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
    gene_half_saturation = gene_half_saturation
  )
  if (!identical(
    colnames(layer1$reaction_expression),
    as.character(layer1$unit_meta$pool_id)
  )) {
    stop(
      "Layer 1 reaction expression and unit metadata are not ordered identically.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}

#' Run directional COMPASS-like LP scoring as a standalone step
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param gem A prepared GEM.
#' @param medium_scenarios Shared medium table.
#' @param outdir Output directory for this stage.
#' @param model_mode Shared union meta-module GEM or shared full GEM.
#' @param layer2_args Additional arguments passed to [rc_run_microcompass()].
#' @param parallel Whether to parallelize shared-model-by-metacell tasks.
#' @param BPPARAM BiocParallel parameter object or `FALSE`.
#' @return The microCOMPASS Layer 2 result.
#' @export
rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(),
    parallel = TRUE,
    BPPARAM = FALSE) {
  model_mode <- match.arg(model_mode)
  if (!inherits(meta_modules, "regcompass_meta_module_step")) {
    stop(
      "`meta_modules` must be the output of `rc_regcompass_step_meta_modules()`.",
      call. = FALSE
    )
  }
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(outdir, "model_cache", model_mode)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  if (!is.list(layer2_args$model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  layer2_args$model_params$cache_dir <- cache_dir
  reserved <- intersect(
    names(layer2_args),
    c(
      "layer1", "gem", "mode", "unit", "reaction_membership",
      "core_reactions", "target_reactions", "medium_scenarios",
      "sample_col", "condition_col", "celltype_col", "BPPARAM", "parallel",
      "penalty_weights"
    )
  )
  if (length(reserved)) {
    stop(
      "`layer2_args` cannot override stepwise workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  global_modules <- meta_modules$global_modules
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = global_modules$global_core_reactions$reaction_id,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      global_modules$global_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      global_modules$global_core_reactions
    } else {
      NULL
    },
    unit = "metacell",
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  answer <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(warning) {
      if (grepl(
        "Metacell-level scores are descriptive pseudo-observations",
        conditionMessage(warning),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!identical(
    colnames(answer$penalty),
    colnames(layer1$reaction_expression)
  )) {
    stop(
      "microCOMPASS and Layer 1 contain different or reordered units.",
      call. = FALSE
    )
  }
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}

#' Assemble rankings and the canonical result after stepwise execution
#'
#' @param metacells Output from [rc_regcompass_step_metacells()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param layer2 Output from [rc_regcompass_step_layer2()].
#' @param gem The prepared GEM used by all stages.
#' @param outdir Final output directory.
#' @param species Species identifier or `"auto"`.
#' @return A canonical RegCompass result list.
#' @export
rc_regcompass_step_results <- function(
    metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse")) {
  if (!inherits(metacells, "regcompass_metacell_step") ||
      !inherits(meta_modules, "regcompass_meta_module_step")) {
    stop(
      "Stepwise results require outputs from the metacell and meta-module stages.",
      call. = FALSE
    )
  }
  species <- .rc_infer_gem_species(gem, species)
  params <- metacells$params
  comparison <- .rc_condition_penalty_comparison(
    layer2,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  conditions <- unique(trimws(as.character(
    metacells$pooled$metacell_meta[[params$condition_col]]
  )))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  result <- list(
    schema_version = "regcompass_v1.7.0_condition_pooled",
    version = "1.7.0",
    species = species,
    model_mode = layer2$model_mode,
    analysis_mode = comparison$analysis_mode,
    pooling_scope = "condition_x_celltype_across_samples",
    input_design = metacells$pooled$input_design,
    metacells = metacells$pooled,
    layer1 = layer1,
    grn_meta_modules = meta_modules$global_modules,
    microcompass = layer2,
    reaction_ranking = comparison$ranking,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    inference_policy = comparison$inference_policy,
    params = list(
      shared_gem = TRUE,
      shared_medium = TRUE,
      n_conditions = length(conditions),
      condition_mode = comparison$analysis_mode,
      metacell_grouping = c(params$condition_col, params$celltype_col),
      samples_mixed_within_condition = any(
        metacells$pooled$metacell_meta$samples_mixed_within_condition
      ),
      sample_weighting = metacells$pooled$sample_weighting,
      biological_sample_minimum = "none; one sample per condition is allowed",
      pando_grouping = c(params$condition_col, params$celltype_col),
      inference_unit = "condition_pooled_metacell_descriptive_only",
      regulatory_alpha = layer1$capacity_params$regulatory_alpha,
      regulatory_state = "ATAC_accessibility_only",
      pando_parameter_source = "RNA_plus_ATAC_condition_x_celltype_fit",
      gpr_promiscuity_mode = "none",
      gpr_and_method = "boltzmann",
      gpr_tau = layer1$capacity_params$tau,
      gpr_or_method = "sum",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      feasibility_completion = "local_fastcore_only",
      penalty_formula = "1/(1+log2(1+E_multiome))",
      reaction_ranking_formula = comparison$ranking_formula,
      execution_mode = "stepwise"
    )
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}
