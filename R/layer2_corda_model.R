# RegCompass integration for the original MATLAB CORDA2 algorithm.

.rc_corda_core_closure_core <- function(
    parent, final, core, target_direction, solver, time_limit,
    flux_threshold) {
  requested <- rc_prepare_directional_targets(
    parent, core, target_direction = target_direction
  )
  parent_diagnostics <- .rc_directional_feasibility(
    parent, requested,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold
  )
  final_diagnostics <- .rc_directional_feasibility(
    final, requested,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold
  )
  names(final_diagnostics)[names(final_diagnostics) == "feasible"] <-
    "final_feasible"
  names(final_diagnostics)[names(final_diagnostics) == "vmax"] <-
    "final_vmax"
  names(final_diagnostics)[names(final_diagnostics) == "solver_status"] <-
    "final_solver_status"
  diagnostics <- merge(
    parent_diagnostics, final_diagnostics,
    by = c("reaction_id", "target_direction"),
    all.x = TRUE, sort = FALSE
  )
  diagnostics$completion_status <- ifelse(
    !diagnostics$feasible,
    "parent_blocked",
    ifelse(
      diagnostics$final_feasible %in% TRUE,
      "corda2_retained",
      "corda2_removed"
    )
  )
  feasible_targets <- diagnostics[
    diagnostics$final_feasible %in% TRUE,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  list(
    requested = requested,
    diagnostics = diagnostics,
    feasible_targets = feasible_targets,
    failed = diagnostics$feasible %in% TRUE &
      !(diagnostics$final_feasible %in% TRUE)
  )
}

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

  included_variables <- reconstruction$included_directional_variables
  if (!length(included_variables)) {
    stop("CORDA2 reconstruction retained no directional reactions.",
         call. = FALSE)
  }
  final <- .rc_corda2_apply_direction_bounds(
    parent, included_variables, split
  )
  closure <- .rc_corda_core_closure(
    parent = parent,
    final = final,
    core = core,
    target_direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = corda_options$flux_threshold
  )

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
  meta$corda2_final_status <- unname(
    reconstruction$final_reaction_status[reaction_id]
  )
  meta$corda2_inclusion_stage <- unname(
    reconstruction$inclusion_stage[reaction_id]
  )
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
  final$target_directions <- closure$feasible_targets
  final$closure_diagnostics <- closure$diagnostics
  final$corda2_reaction_evidence <- reaction_evidence
  final$corda_reaction_evidence <- reaction_evidence
  final$corda_task_diagnostics <- reconstruction$task_diagnostics
  final$corda_execution <- list(original_matlab = list(
    solver_runtime = if (isTRUE(reconstruction$solver_performance$persistent_solver)) {
      "highs_persistent_cpp"
    } else {
      "one_shot"
    },
    target_parallelism = FALSE
  ))
  final$corda_stage1_HCtoMC <- reconstruction$HCtoMC
  final$corda_stage1_HCtoNC <- reconstruction$HCtoNC
  final$corda_stage2_MCtoNC <- reconstruction$MCtoNC
  final$corda_rescue <- reconstruction$rescue
  final$corda_reconstruction <- reconstruction
  final$target_status <- if (any(closure$failed)) {
    "core_direction_removed_by_corda2"
  } else if (!nrow(closure$feasible_targets)) {
    "no_feasible_core_direction"
  } else {
    "ok"
  }
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
    reconstruction_direction_policy =
      "original_CORDA2_opposite_direction_closed_and_directional_merge",
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
    source_semantics = reconstruction$source_semantics
  )

  final <- .rc_corda2_apply_target_flux(
    model = final,
    flux_threshold = corda_options$flux_threshold,
    strict = strict,
    cell_type = cell_type
  )
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

# Progress-aware entry point; the algorithm remains in the core above.
.rc_corda_core_closure <- function(...) {
  task <- .rc_layer2_progress_state$current_task
  if (!is.null(task) && identical(task$route, "corda2")) {
    .rc_layer2_current_task_event(
      "corda2_target_closure", 8L,
      "validating retained core directions in the reconstructed GEM"
    )
  }
  do.call(.rc_corda_core_closure_core, list(...))
}
