# RegCompass integration for the original MATLAB CORDA2 algorithm.

.rc_complete_celltype_medium_corda_gem_core <- function(
    gem, reaction_membership, core_reactions, cell_type,
    reaction_evidence, corda_options, medium_table = NULL,
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300, fastcore_epsilon = 1e-4,
    max_support_reactions = 2000, strict = TRUE) {
  target_direction <- match.arg(target_direction)
  forbidden_roles <- c("demand", "sink", "artificial_support")
  parent <- .rc_corda_parent(
    gem = gem,
    medium_table = medium_table,
    condition = NULL,
    forbidden_roles = forbidden_roles,
    solver = solver,
    time_limit = time_limit
  )
  validated <- rc_validate_gem(parent)
  module_reactions <- intersect(
    unique(as.character(reaction_membership$reaction_id)),
    validated$reactions
  )
  core <- intersect(
    unique(as.character(core_reactions$reaction_id)), module_reactions
  )
  if (!length(module_reactions) || !length(core)) {
    stop("CORDA2 requires valid module and core reactions.", call. = FALSE)
  }
  classes <- .rc_corda_classify_reactions(
    parent_reactions = validated$reactions,
    module_reactions = module_reactions,
    core_reactions = core,
    reaction_evidence = reaction_evidence,
    medium_confidence_threshold =
      corda_options$medium_confidence_threshold,
    negative_confidence_threshold =
      corda_options$negative_confidence_threshold,
    include_evidence_outside_modules =
      corda_options$include_evidence_outside_modules,
    max_medium_confidence_reactions =
      corda_options$max_medium_confidence_reactions
  )
  split <- .rc_corda2_split_original(parent)
  reconstruction <- .rc_corda_build_three_stage(
    split = split,
    classes = classes,
    options = corda_options,
    solver = solver,
    time_limit = time_limit
  )
  reconstruction$source_fidelity <- "original_MATLAB_CORDA2"
  reconstruction$solver_time_limit <- time_limit
  stage_parallel <- identical(
    reconstruction$parallel_execution_policy,
    "stage_barrier_parallel_targets_deterministic_ordered_reduce"
  )
  solver_state_scope <- paste(
    "fresh solver engine per directional target; persistent simplex basis",
    "reuse only inside that target's maximize/dependency iterations"
  )
  reconstruction$solver_state_scope <- solver_state_scope
  if (!stage_parallel) {
    reconstruction$parallel_execution_policy <-
      "serial_original_target_order_target_isolated_solver_state"
  }

  included_variables <- unique(as.character(
    reconstruction$included_directional_variables %||% character()
  ))
  final <- .rc_corda2_apply_direction_bounds(
    parent = parent,
    included_variables = included_variables,
    split = split,
    core_reactions = core
  )
  missing_core <- setdiff(core, colnames(final$S))
  if (length(missing_core)) {
    stop(
      "CORDA2 finalization failed to retain required core reactions: ",
      paste(utils::head(missing_core, 10L), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  selected_by_corda2 <- unique(as.character(
    split$variable_to_reaction[
      intersect(included_variables, split$variable_order)
    ]
  ))
  selected_by_corda2 <- selected_by_corda2[
    !is.na(selected_by_corda2) & nzchar(selected_by_corda2)
  ]
  forced_core <- setdiff(core, selected_by_corda2)

  meta <- final$reaction_meta
  if (is.null(meta)) {
    meta <- data.frame(
      reaction_id = colnames(final$S), stringsAsFactors = FALSE
    )
  }
  reaction_id <- as.character(meta$reaction_id)
  evidence <- reaction_evidence[
    match(reaction_id, as.character(reaction_evidence$reaction_id)),
    , drop = FALSE
  ]
  meta$merged_meta_module_member <- reaction_id %in% module_reactions
  meta$core_reaction <- reaction_id %in% core
  meta$core_structural_retention <- ifelse(
    meta$core_reaction, "required", "not_core"
  )
  meta$celltype_fastcore_support <- FALSE
  meta$celltype_corda2_included <- TRUE
  meta$celltype_corda_included <- TRUE
  meta$support_only <- !reaction_id %in% module_reactions
  meta$corda2_initial_confidence <- unname(
    reconstruction$initial_reaction_confidence[reaction_id]
  )
  meta$corda2_initial_class <- unname(
    classes$initial_confidence[reaction_id]
  )
  meta$corda2_final_confidence <- unname(
    reconstruction$final_reaction_confidence[reaction_id]
  )
  directional_status <- unname(
    reconstruction$final_reaction_status[reaction_id]
  )
  meta$corda2_directional_status <- directional_status
  meta$corda2_final_status <- ifelse(
    meta$core_reaction & !(directional_status %in% "included"),
    "core_forced_retained",
    directional_status
  )
  meta$corda2_inclusion_stage <- unname(
    reconstruction$inclusion_stage[reaction_id]
  )
  meta$corda2_inclusion_stage[
    meta$core_reaction & is.na(meta$corda2_inclusion_stage)
  ] <- "core_structural_backbone"
  meta$corda2_evidence_score <- unname(classes$evidence_score[reaction_id])
  meta$corda2_rna_percentile <- suppressWarnings(as.numeric(
    evidence$rna_percentile
  ))
  meta$corda2_multiome_percentile <- suppressWarnings(as.numeric(
    evidence$multiome_percentile
  ))
  meta$corda2_regulatory_support <- suppressWarnings(as.numeric(
    evidence$regulatory_support
  ))
  meta$corda_initial_confidence <- meta$corda2_initial_confidence
  meta$corda_final_confidence <- meta$corda2_final_confidence
  meta$corda_inclusion_stage <- meta$corda2_inclusion_stage
  meta$corda_evidence_score <- meta$corda2_evidence_score
  final$reaction_meta <- meta
  final$sample_id <- cell_type
  final$grn_module_id <- paste0("CELLTYPE_MEDIUM_UNION_GEM::", cell_type)
  final$cell_type <- cell_type
  final$required_core_reactions <- core
  final$corda2_reaction_evidence <- reaction_evidence
  final$corda_reaction_evidence <- reaction_evidence
  final$corda_task_diagnostics <- reconstruction$task_diagnostics
  final$corda_execution <- list(original_matlab = list(
    solver_runtime = if (isTRUE(reconstruction$solver_performance$persistent_solver)) {
      "highs_persistent_cpp"
    } else {
      "one_shot"
    },
    solver_state_scope = solver_state_scope,
    target_parallelism = stage_parallel,
    parallel_scope = if (stage_parallel) {
      "directional_targets_within_each_original_corda2_stage"
    } else {
      "serial_original_target_order_with_target_isolated_solver_state"
    },
    stage_barrier = stage_parallel,
    ordered_reduce = stage_parallel,
    worker_lifecycle = if (stage_parallel) {
      "fresh_pool_each_stage_stop_full_gc_before_next_stage"
    } else {
      "no_worker_pool_fresh_solver_engine_per_target"
    },
    closure_parallelism = "none_post_reconstruction"
  ))
  final$corda_stage1_HCtoMC <- reconstruction$HCtoMC
  final$corda_stage1_HCtoNC <- reconstruction$HCtoNC
  final$corda_stage2_MCtoNC <- reconstruction$MCtoNC
  final$corda_rescue <- reconstruction$rescue
  final$corda_reconstruction <- reconstruction
  final$is_union_gem <- TRUE
  final$union_gem_scope <-
    "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"

  original_args <- list(
    MCxNCthresh = corda_options$MCxNCthresh,
    constraint = corda_options$constraint,
    constrainby = corda_options$constrainby,
    om = corda_options$om,
    ci = corda_options$ci
  )
  initial <- classes$initial_confidence
  included <- colnames(final$S)
  final$build_params <- list(
    strategy = "celltype_medium_original_matlab_corda2",
    cell_type = cell_type,
    algorithm = reconstruction$algorithm,
    completion_stage = "original_CORDA2_after_confidence_mapping",
    evidence_schema = "regcompass_corda_reaction_evidence_v2",
    confidence_mapping = classes$confidence_contract,
    n_celltype_biological_reactions = length(included),
    n_celltype_fastcore_support_reactions = 0L,
    n_high_confidence_reactions = sum(initial == "HC"),
    n_module_medium_confidence_reactions = sum(initial == "MC_module"),
    n_evidence_medium_confidence_reactions = sum(initial == "MC_evidence"),
    n_negative_confidence_reactions = sum(initial == "NC"),
    n_other_reactions = sum(initial == "OT"),
    n_corda_included_reactions = length(included),
    n_corda2_directionally_selected_reactions = length(selected_by_corda2),
    n_core_reactions = length(core),
    n_core_forced_retained = length(forced_core),
    n_corda_included_initial_HC = sum(
      included %in% names(initial)[initial == "HC"]
    ),
    n_corda_included_initial_MC = sum(
      included %in% names(initial)[grepl("^MC", initial)]
    ),
    n_corda_included_initial_NC = sum(
      included %in% names(initial)[initial == "NC"]
    ),
    n_corda_included_initial_OT = sum(
      included %in% names(initial)[initial == "OT"]
    ),
    n_stage1_associated = length(reconstruction$stage1_associated),
    n_stage2_promoted_nc = length(reconstruction$stage2_promoted_nc),
    n_stage2_promoted_mc = length(reconstruction$stage2_promoted_mc),
    n_stage3_associated_ot = length(reconstruction$stage3_associated_ot),
    scoring_target_direction = target_direction,
    reconstruction_direction_policy = paste(
      "original CORDA2 directional state machine internally; retain a",
      "reaction if either split direction is selected; retain every core",
      "reaction unconditionally and restore medium-constrained parent bounds"
    ),
    core_retention_policy = "immutable_structural_backbone",
    closure_policy = "none_post_reconstruction_compass_vmax_only",
    corda2_args = original_args,
    corda2_solver_time_limit = time_limit,
    fastcore_epsilon_used = FALSE,
    max_support_reactions_ignored_for_corda2 = max_support_reactions,
    strict_requested = strict,
    strict_used_for_reconstruction = FALSE,
    corda_options = corda_options,
    included_reactions = included,
    included_directional_variables = included_variables,
    stage_update_policy = reconstruction$stage_update_policy,
    parallel_execution_policy = reconstruction$parallel_execution_policy,
    solver_state_scope = solver_state_scope,
    source_semantics = reconstruction$source_semantics
  )
  final$corda2_contract <- list(
    implementation = "original MATLAB CORDA2.m semantics",
    reference_repository = corda_options$reference_repository,
    reference_file = corda_options$reference_file,
    adjustable_args = original_args,
    fixed_internal = c(
      fluxThreshold = corda_options$flux_threshold,
      baselineCost = corda_options$baseline_cost,
      outputBound = corda_options$output_bound
    ),
    solver_time_limit = time_limit,
    stage_update_policy = reconstruction$stage_update_policy,
    parallel_execution_policy = reconstruction$parallel_execution_policy,
    solver_state_scope = solver_state_scope,
    source_semantics = reconstruction$source_semantics,
    post_reconstruction_direction_merge = paste(
      "retain reaction if either direction is selected; force all core",
      "reactions; restore medium-constrained parent bounds"
    ),
    core_retention = "all_core_reactions_structurally_required",
    post_reconstruction_closure =
      "none_microcompass_scoring_computes_directional_vmax_once"
  )

  final <- .rc_corda2_prepare_scoring_targets(
    model = final,
    core_reactions = core,
    target_direction = target_direction,
    strict = strict,
    cell_type = cell_type
  )
  progress_state <- get0(
    ".rc_layer2_progress_state", mode = "environment", inherits = TRUE
  )
  task <- if (is.environment(progress_state)) {
    progress_state$current_task
  } else {
    NULL
  }
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "corda2_scoring_targets_ready", 8L,
      detail = paste0(
        "core_reactions=", length(core),
        "; directional_targets=", nrow(final$target_directions),
        "; post_build_lp_solves=0"
      )
    )
  }

  final <- .rc_corda_attach_parent_contract(
    model = final,
    parent = parent,
    fastcore_epsilon = fastcore_epsilon,
    forbidden_roles = forbidden_roles
  )
  .rc_finalize_corda_union_model(final, cell_type = cell_type)
}

# Progress-aware entry point; the algorithm remains in the core above.
.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  context <- .rc_layer2_task_context(
    cell_type = args$cell_type %||% args[[4L]],
    medium_scenario = .rc_layer2_medium_id(args$medium_table),
    route = "corda2"
  )
  previous <- .rc_layer2_task_push(context, "corda2", 9L)
  on.exit(.rc_layer2_task_pop(previous), add = TRUE)
  tryCatch({
    answer <- do.call(
      .rc_complete_celltype_medium_corda_gem_core,
      args
    )
    .rc_layer2_current_task_event(
      "corda2_model_ready", 9L,
      detail = paste0(
        "included_reactions=", ncol(answer$S),
        "; lp_solves=",
        answer$corda_reconstruction$solver_performance$n_solves %||% NA_integer_
      ),
      status = "complete"
    )
    answer
  }, error = function(error) {
    .rc_layer2_current_task_event(
      "corda2_task_error",
      .rc_layer2_task_last_step(context),
      conditionMessage(error), status = "error"
    )
    stop(error)
  })
}
