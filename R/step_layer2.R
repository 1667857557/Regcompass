# Layer 2 shared-structure validation and comparison helpers.
.rc_layer2_direction_contract <- function(x) {
  tab <- as.data.frame(x$lp_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "medium_scenario", "row_id"
  )
  if (!all(required %in% colnames(tab))) {
    stop("Layer 2 direction diagnostics are incomplete.", call. = FALSE)
  }
  tab <- unique(tab[, required, drop = FALSE])
  tab[order(
    tab$medium_scenario, tab$reaction_id, tab$target_direction, tab$row_id
  ), , drop = FALSE]
}

.rc_assert_layer2_shared_contract <- function(primary, candidate, label) {
  matrices <- c("penalty", "vmax", "feasible", "evaluated")
  for (name in matrices) {
    if (!identical(dimnames(primary[[name]]), dimnames(candidate[[name]]))) {
      stop("Layer 2 ", label, " changed ", name, " ordering.",
           call. = FALSE)
    }
  }
  if (!identical(
        primary$structural_model_contract,
        candidate$structural_model_contract
      ) ||
      !identical(
        .rc_layer2_direction_contract(primary),
        .rc_layer2_direction_contract(candidate)
      ) ||
      !identical(primary$medium_scenarios, candidate$medium_scenarios) ||
      !identical(primary$unit_meta, candidate$unit_meta)) {
    stop(
      "Layer 2 ", label,
      " did not reuse the exact GEM, bounds, media, directions, and units.",
      call. = FALSE
    )
  }
  finite <- is.finite(primary$vmax) & is.finite(candidate$vmax)
  if (any(is.finite(primary$vmax) != is.finite(candidate$vmax)) ||
      any(abs(primary$vmax[finite] - candidate$vmax[finite]) >
          1e-10 * pmax(1, abs(primary$vmax[finite])))) {
    stop("Layer 2 ", label, " changed structural target vmax.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_layer2_comparison_table <- function(
    layer2, layer1, condition_col, celltype_col) {
  penalty <- layer2$penalty_common
  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  unit_meta <- layer2$unit_meta
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else {
    as.character(unit_meta$pool_id)
  }
  unit_meta <- unit_meta[match(colnames(penalty), unit_id), , drop = FALSE]
  grid <- expand.grid(
    row_index = seq_len(nrow(penalty)),
    unit_index = seq_len(ncol(penalty)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  reaction <- row_meta$reaction_id[grid$row_index]
  unit <- colnames(penalty)[grid$unit_index]
  index_matrix <- function(x) {
    as.numeric(x[cbind(grid$row_index, grid$unit_index)])
  }
  reaction_unit_matrix <- function(x) {
    if (!is.numeric(x) || is.null(dim(x)) ||
        !all(unique(row_meta$reaction_id) %in% rownames(x)) ||
        !identical(colnames(x), colnames(penalty))) {
      stop("Layer 1 reaction diagnostics are not aligned.", call. = FALSE)
    }
    as.numeric(x[cbind(
      match(reaction, rownames(x)),
      grid$unit_index
    )])
  }
  lookup_flag <- function(table, flag) {
    if (!is.data.frame(table) ||
        !all(c("reaction_id", "cell_type", flag) %in% colnames(table))) {
      return(rep(NA, nrow(grid)))
    }
    key <- paste(
      as.character(table$reaction_id),
      as.character(table$cell_type),
      sep = "\001"
    )
    requested <- paste(
      reaction,
      as.character(unit_meta[[celltype_col]][grid$unit_index]),
      sep = "\001"
    )
    as.logical(table[[flag]][match(requested, key)])
  }
  omega <- layer2$params$omega
  vmax <- index_matrix(layer2$vmax)
  common <- index_matrix(layer2$penalty_common)
  normalized <- common / (omega * vmax)
  normalized[!is.finite(common) | !is.finite(vmax) | vmax <= 0] <- NA_real_
  common_fraction <- reaction_unit_matrix(
    layer1$reaction_common_support_fraction
  )
  full_fraction <- reaction_unit_matrix(
    layer1$reaction_condition_full_support_fraction
  )
  data.frame(
    reaction_id = reaction,
    direction = row_meta$target_direction[grid$row_index],
    medium = row_meta$medium_scenario[grid$row_index],
    cell_type =
      as.character(unit_meta[[celltype_col]][grid$unit_index]),
    condition =
      as.character(unit_meta[[condition_col]][grid$unit_index]),
    metacell_id = unit,
    penalty_rna_only = index_matrix(layer2$penalty_rna_only),
    penalty_common_oof = common,
    penalty_condition_full_oof =
      index_matrix(layer2$penalty_condition_full),
    penalty_unique_increment =
      index_matrix(layer2$penalty_unique_increment),
    penalty_per_target_flux = normalized,
    vmax = vmax,
    projection_oof_available =
      is.finite(common) & is.finite(common_fraction) &
      common_fraction >= 1 - sqrt(.Machine$double.eps),
    common_support_fraction = common_fraction,
    condition_full_support_fraction = full_fraction,
    depth_sensitivity_flag = lookup_flag(
      layer1$depth_diagnostics$reaction_depth_sensitivity,
      "depth_sensitivity_flag"
    ),
    zero_support_sensitive = lookup_flag(
      layer1$reaction_zero_support_sensitivity,
      "zero_support_sensitive"
    ),
    link_saturation_sensitive = lookup_flag(
      layer1$reaction_link_saturation_sensitivity,
      "link_saturation_sensitive"
    ),
    alpha = layer1$capacity_params$regulatory_alpha,
    inference_class = "metacell_statistical_unit_within_dataset",
    comparability_class = "common_support_primary",
    stringsAsFactors = FALSE
  )
}

#' Build medium-specific structural models and run directional LP scoring
#'
#' With `model_mode = "meta_module_gem"`, this stage is the only place where
#' FASTCORE is applied. For each medium scenario it constructs one union GEM
#' from the merged biological meta-module catalogue plus global FASTCORE
#' support, then reuses that exact model for every condition, evidence route,
#' and metacell. Only union-GEM construction accepts
#' `model_params$completion_time_limit`; scoring LPs have no time-limit
#' parameter.
#'
#' @export
rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer2", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  model_mode <- match.arg(model_mode)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c(
    "model_params", "omega", "target_direction", "solver", "flux_threshold"
  )
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop(
      "Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
      ". Scoring `time_limit` has been removed. Use only ",
      "`layer2_args$model_params$completion_time_limit` to limit global ",
      "FASTCORE union-GEM construction.",
      call. = FALSE
    )
  }
  params <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  if (!is.list(layer2_args$model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  allowed_model_params <- c(
    "completion_time_limit", "fastcore_epsilon",
    "max_support_reactions", "strict"
  )
  unknown_model_params <- setdiff(
    names(layer2_args$model_params), allowed_model_params
  )
  if (length(unknown_model_params)) {
    stop(
      "Unsupported `layer2_args$model_params`: ",
      paste(unknown_model_params, collapse = ", "),
      ". Only union-GEM construction controls are accepted; ",
      "`time_limit` is not a scoring or model parameter.",
      call. = FALSE
    )
  }
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
    stop(
      "`layer2_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  .rc_require_lp_solver(solver)
  catalogue <- meta_modules$merged_modules
  if (!is.list(catalogue) ||
      !is.data.frame(catalogue$merged_core_reactions) ||
      !is.data.frame(catalogue$merged_reaction_membership)) {
    stop(
      "The merged biological meta-module catalogue is incomplete.",
      call. = FALSE
    )
  }
  targets <- unique(as.character(
    catalogue$merged_core_reactions$reaction_id
  ))
  missing_expression <- setdiff(
    targets, rownames(layer1$reaction_expression)
  )
  if (length(missing_expression)) {
    stop(
      "Merged core reactions are absent from Layer 1 expression: ",
      paste(utils::head(missing_expression, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = targets,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_core_reactions
    } else {
      NULL
    },
    unit = "metacell",
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  answer <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(w) {
      if (grepl(
        "Metacells are valid within-dataset statistical units",
        conditionMessage(w), fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  run_comparator <- function(expression_field, label) {
    expression <- if (is.character(expression_field)) {
      layer1[[expression_field]]
    } else {
      expression_field
    }
    if (!is.numeric(expression) || is.null(dim(expression)) ||
        !identical(
          dimnames(expression),
          dimnames(layer1$reaction_expression)
        )) {
      stop(
        "Layer 1 comparator `", label,
        "` is missing or is not aligned with the common OOF expression.",
        call. = FALSE
      )
    }
    comparator_layer1 <- layer1
    comparator_layer1$reaction_expression <- expression
    comparator_args <- c(defaults, layer2_args)
    comparator_args$layer1 <- comparator_layer1
    comparator_args$model_cache_override <- answer$shared_model_cache
    result <- withCallingHandlers(
      do.call(.rc_run_microcompass_engine, comparator_args),
      warning = function(w) {
        if (grepl(
          "Metacells are valid within-dataset statistical units",
          conditionMessage(w), fixed = TRUE
        )) invokeRestart("muffleWarning")
      }
    )
    .rc_assert_layer2_shared_contract(answer, result, label)
    result
  }
  condition_full <- run_comparator(
    "reaction_expression_condition_full_oof",
    "condition-full OOF comparator"
  )
  rna_only <- run_comparator(
    "reaction_expression_rna_only",
    "RNA-only comparator"
  )
  depth_matched <- run_comparator(
    "reaction_expression_depth_matched_rna",
    "depth-matched RNA comparator"
  )
  common_depth_interval <- run_comparator(
    "reaction_expression_common_depth_interval_rna",
    "common-depth-interval RNA comparator"
  )
  alpha_paths <- lapply(
    names(layer1$reaction_expression_alpha_sensitivity),
    function(alpha_name) {
      run_comparator(
        layer1$reaction_expression_alpha_sensitivity[[alpha_name]],
        paste0("alpha=", alpha_name, " comparator")
      )
    }
  )
  names(alpha_paths) <- names(
    layer1$reaction_expression_alpha_sensitivity
  )
  answer$penalty_common <- answer$penalty
  answer$penalty_condition_full <- condition_full$penalty
  answer$penalty_rna_only <- rna_only$penalty
  answer$penalty_depth_matched_rna <- depth_matched$penalty
  answer$penalty_common_depth_interval_rna <-
    common_depth_interval$penalty
  answer$penalty_alpha_sensitivity <- lapply(
    alpha_paths, `[[`, "penalty"
  )
  answer$penalty_unique_increment <-
    answer$penalty_condition_full - answer$penalty_common
  answer$penalty_grn_total_increment <-
    answer$penalty_condition_full - answer$penalty_rna_only
  answer$score_common_display_only <- answer$score
  answer$score_condition_full_display_only <- condition_full$score
  answer$score_rna_only_display_only <- rna_only$score
  answer$score_depth_matched_rna_display_only <- depth_matched$score
  answer$score_common_depth_interval_rna_display_only <-
    common_depth_interval$score
  answer$score_alpha_sensitivity_display_only <- lapply(
    alpha_paths, `[[`, "score"
  )
  answer$comparison_contract <- list(
    primary = "penalty_common",
    exploratory = "penalty_condition_full",
    control = "penalty_rna_only",
    depth_sensitivity = c(
      "penalty_depth_matched_rna",
      "penalty_common_depth_interval_rna"
    ),
    alpha_sensitivity = paste0(
      "penalty_alpha_sensitivity[[", names(alpha_paths), "]]"
    ),
    unique_increment =
      "penalty_condition_full - penalty_common",
    grn_total_increment =
      "penalty_condition_full - penalty_rna_only",
    recomputation =
      "complete GPR, penalty, and LP rerun for each evidence route",
    exact_shared_structure = TRUE,
    structural_model_contract = answer$structural_model_contract,
    effect_size_basis = "penalty / (omega * vmax)",
    ecdf_effect_size_eligible = FALSE
  )
  condition_full$shared_model_cache <- NULL
  rna_only$shared_model_cache <- NULL
  depth_matched$shared_model_cache <- NULL
  common_depth_interval$shared_model_cache <- NULL
  alpha_paths <- lapply(alpha_paths, function(x) {
    x$shared_model_cache <- NULL
    x
  })
  answer$comparison_paths <- list(
    condition_full_oof = condition_full,
    rna_only = rna_only,
    depth_matched_rna = depth_matched,
    common_depth_interval_rna = common_depth_interval,
    alpha_sensitivity = alpha_paths
  )
  answer$comparison_table <- .rc_layer2_comparison_table(
    answer, layer1, params$condition_col, params$celltype_col
  )
  answer$workflow_params <- params
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  answer$source_core_reactions <- catalogue$merged_core_reactions
  answer$source_merged_reaction_membership <-
    catalogue$merged_reaction_membership
  answer$union_gem_policy <- if (identical(model_mode, "meta_module_gem")) {
    "one medium-specific union GEM; single global FASTCORE completion"
  } else {
    "shared full GEM; no union-GEM reconstruction"
  }
  class(answer) <- c("regcompass_layer2_step", "list")
  .rc_validate_layer2_stage(
    answer,
    layer1 = layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer2"
  )
  answer <- .rc_step_monitor_finish(answer, monitor)
  rc_export_microcompass(answer, outdir)
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}
