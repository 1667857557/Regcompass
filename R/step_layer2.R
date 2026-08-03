# Layer 2 shared-structure validation and condition-full scoring.

.rc_layer2_direction_contract <- function(x) {
  tab <- as.data.frame(x$lp_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "medium_scenario", "row_id"
  )
  if (identical(as.character(x$model_mode), "meta_module_gem")) {
    required <- c("cell_type", required)
  }
  if (!all(required %in% colnames(tab))) {
    stop("Layer 2 direction diagnostics are incomplete.", call. = FALSE)
  }
  tab <- unique(tab[, required, drop = FALSE])
  order_cols <- intersect(
    c("cell_type", "medium_scenario", "reaction_id",
      "target_direction", "row_id"),
    colnames(tab)
  )
  tab[do.call(order, tab[order_cols]), , drop = FALSE]
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
  penalty <- layer2$penalty
  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  unit_meta <- layer2$unit_meta
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else {
    as.character(unit_meta$pool_id)
  }
  unit_meta <- unit_meta[match(colnames(penalty), unit_id), , drop = FALSE]
  unit_celltype <- as.character(unit_meta[[celltype_col]])
  grid <- expand.grid(
    row_index = seq_len(nrow(penalty)),
    unit_index = seq_len(ncol(penalty)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scoped_celltype <- row_meta$cell_type[grid$row_index]
  matching_scope <- is.na(scoped_celltype) |
    scoped_celltype == unit_celltype[grid$unit_index]
  grid <- grid[matching_scope, , drop = FALSE]
  if (!nrow(grid)) {
    stop("No Layer 2 rows match their cell-type units.", call. = FALSE)
  }
  reaction <- row_meta$reaction_id[grid$row_index]
  unit <- colnames(penalty)[grid$unit_index]
  index_matrix <- function(x) {
    as.numeric(x[cbind(grid$row_index, grid$unit_index)])
  }
  reaction_unit_matrix <- function(x) {
    if (!is.numeric(x) || is.null(dim(x)) ||
        !all(unique(reaction) %in% rownames(x)) ||
        !identical(colnames(x), colnames(penalty))) {
      stop("Layer 1 reaction diagnostics are not aligned.", call. = FALSE)
    }
    as.numeric(x[cbind(match(reaction, rownames(x)), grid$unit_index)])
  }
  omega <- layer2$params$omega
  vmax <- index_matrix(layer2$vmax)
  primary <- index_matrix(layer2$penalty_condition_full_oof)
  normalized <- primary / (omega * vmax)
  normalized[!is.finite(primary) | !is.finite(vmax) | vmax <= 0] <- NA_real_
  data.frame(
    row_id = rownames(penalty)[grid$row_index],
    reaction_id = reaction,
    direction = row_meta$target_direction[grid$row_index],
    medium = row_meta$medium_scenario[grid$row_index],
    cell_type = unit_celltype[grid$unit_index],
    condition = as.character(
      unit_meta[[condition_col]][grid$unit_index]
    ),
    metacell_id = unit,
    penalty_condition_full_oof = primary,
    penalty_common_oof = index_matrix(layer2$penalty_common_oof),
    penalty_condition_unique_increment = index_matrix(
      layer2$penalty_condition_unique_increment
    ),
    penalty_rna_only = index_matrix(layer2$penalty_rna_only),
    penalty_per_target_flux = normalized,
    vmax = vmax,
    condition_full_oof_available = is.finite(primary),
    condition_full_support_fraction = reaction_unit_matrix(
      layer1$reaction_condition_full_support_fraction
    ),
    common_support_fraction = reaction_unit_matrix(
      layer1$reaction_common_support_fraction
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    comparability_class =
      "same_celltype_conditions_on_one_celltype_medium_union_gem",
    stringsAsFactors = FALSE
  )
}

.rc_compact_meta_modules_for_layer2 <- function(meta_modules) {
  embedded <- meta_modules$condition_modules
  if (is.list(embedded) && is.data.frame(embedded$reaction_membership)) {
    meta_modules$condition_modules <- list(
      schema_version = "regcompass_transient_condition_modules_summary_v1",
      embedded = FALSE,
      n_reaction_membership = nrow(embedded$reaction_membership)
    )
    class(meta_modules$condition_modules) <- c(
      "regcompass_external_condition_modules", "list"
    )
  }
  meta_modules
}

.rc_run_microcompass_monitored <- function(..., progress_monitor = NULL) {
  args <- list(...)
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_engine_start",
    "starting union-model construction and directional scoring",
    context = list(parallel = isTRUE(args$parallel %||% TRUE))
  )
  answer <- do.call(rc_run_microcompass, args)
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_engine_complete",
    "completed union-model construction and directional scoring"
  )
  answer
}

.rc_run_microcompass_engine_monitored <- function(
    ..., progress_monitor = NULL) {
  args <- list(...)
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_control_start",
    "starting shared-model control scoring"
  )
  answer <- do.call(.rc_run_microcompass_engine, args)
  .rc_step_monitor_event(
    progress_monitor,
    "layer2_control_complete",
    "completed shared-model control scoring"
  )
  answer
}

#' Build cell-type and medium structural models and score directional LPs
#'
#' With `model_mode = "meta_module_gem"`, conditions are unioned only within the
#' same cell type. One union GEM and one independent FASTCORE completion are
#' created for every cell-type and medium combination. Conditions and metacells
#' share a model only when their cell type matches. Historical `_oof`, `common`,
#' and `condition_unique` fields remain compatibility aliases; RNA-only scoring
#' is an interpretation control. The `full_gem` route uses a separate engine.
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
  meta_modules <- .rc_compact_meta_modules_for_layer2(meta_modules)
  invisible(gc(verbose = FALSE, full = TRUE))
  .rc_require_stage_class(
    meta_modules,
    "regcompass_meta_module_step",
    "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c(
    "model_params", "omega", "target_direction", "solver",
    "flux_threshold"
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
    layer1,
    workflow_params = params,
    gem = gem,
    argument = "layer1"
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
      ". Only union-GEM construction controls are accepted.",
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
    stop("The merged biological meta-module catalogue is incomplete.",
         call. = FALSE)
  }
  required_catalogue_cols <- c(params$celltype_col, "reaction_id")
  if (!is.list(catalogue$cell_type_catalogues) ||
      !all(required_catalogue_cols %in%
           colnames(catalogue$merged_core_reactions)) ||
      !all(required_catalogue_cols %in%
           colnames(catalogue$merged_reaction_membership)) ||
      !identical(catalogue$merge_scope, "cell_type") ||
      isTRUE(catalogue$cross_celltype_merge)) {
    stop("Meta-modules are not partitioned by cell type.", call. = FALSE)
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
    target_reactions = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_core_reactions
    } else {
      targets
    },
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
    do.call(
      .rc_run_microcompass_monitored,
      c(defaults, layer2_args, list(progress_monitor = monitor))
    ),
    warning = function(w) {
      if (grepl(
        "Metacells are valid within-dataset statistical units",
        conditionMessage(w),
        fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  run_control <- function(expression_field, label) {
    expression <- layer1[[expression_field]]
    if (!is.numeric(expression) || is.null(dim(expression)) ||
        !identical(
          dimnames(expression),
          dimnames(layer1$reaction_expression)
        )) {
      stop(
        "Layer 1 control `", label,
        "` is missing or is not aligned with the primary expression.",
        call. = FALSE
      )
    }
    control_layer1 <- layer1
    control_layer1$reaction_expression <- expression
    control_args <- c(defaults, layer2_args)
    control_args$layer1 <- control_layer1
    control_args$model_cache_override <- answer$shared_model_cache
    result <- withCallingHandlers(
      do.call(
        .rc_run_microcompass_engine_monitored,
        c(control_args, list(progress_monitor = monitor))
      ),
      warning = function(w) {
        if (grepl(
          "Metacells are valid within-dataset statistical units",
          conditionMessage(w),
          fixed = TRUE
        )) invokeRestart("muffleWarning")
      }
    )
    .rc_assert_layer2_shared_contract(answer, result, label)
    result
  }
  common <- run_control(
    "reaction_expression_common_oof",
    "common-support compatibility route"
  )
  rna_only <- run_control(
    "reaction_expression_rna_only",
    "RNA-only control"
  )
  answer$schema_version <- "regcompass_regulatory_layer2_v2"
  answer$penalty_condition_full_oof <- answer$penalty
  answer$penalty_common_oof <- common$penalty
  answer$penalty_condition_unique_increment <-
    answer$penalty_condition_full_oof - answer$penalty_common_oof
  answer$penalty_rna_only <- rna_only$penalty
  answer$score_condition_full_oof_display_only <- answer$score
  answer$score_common_oof_display_only <- common$score
  answer$score_rna_only_display_only <- rna_only$score
  answer$comparison_contract <- list(
    primary = "penalty_condition_full_oof",
    common_component = "penalty_common_oof",
    condition_unique_increment =
      "penalty_condition_full_oof - penalty_common_oof",
    rna_control = "penalty_rna_only",
    nonestimable_edge_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    removed_guardrails = c(
      "depth_matching",
      "common_depth_restriction",
      "alpha_sensitivity",
      "zero_support_sensitivity",
      "link_saturation_propagation"
    ),
    exact_shared_structure = TRUE,
    structural_model_contract = answer$structural_model_contract,
    effect_size_basis = "penalty / (omega * vmax)",
    ecdf_effect_size_eligible = FALSE
  )
  common$shared_model_cache <- NULL
  rna_only$shared_model_cache <- NULL
  answer$comparison_paths <- list(
    common_support = common,
    rna_only = rna_only
  )
  answer$comparison_table <- .rc_layer2_comparison_table(
    answer,
    layer1,
    params$condition_col,
    params$celltype_col
  )
  answer$workflow_params <- params
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  answer$source_core_reactions <- catalogue$merged_core_reactions
  answer$source_merged_reaction_membership <-
    catalogue$merged_reaction_membership
  answer$union_gem_policy <- if (identical(model_mode, "meta_module_gem")) {
    paste(
      "one union GEM per cell type and medium; FASTCORE runs",
      "independently within each cell type"
    )
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
