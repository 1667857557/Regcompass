# Build one evidence-maximizing CORDA-like union GEM.

.rc_complete_celltype_medium_corda_like_gem <- function(
    gem, reaction_membership, core_reactions, cell_type,
    reaction_evidence, corda_options, medium_table = NULL,
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300, fastcore_epsilon = 1e-4,
    max_support_reactions = 2000, strict = TRUE) {
  target_direction <- match.arg(target_direction)
  parent <- .rc_fastcore_parent(
    gem,
    medium_table = medium_table,
    condition = NULL,
    solver = solver,
    time_limit = time_limit,
    fastcore_epsilon = fastcore_epsilon
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
    stop("CORDA-like completion requires valid module and core reactions.",
         call. = FALSE)
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
  support_costs <- .rc_corda_support_costs(
    validated$reactions,
    classes,
    other_penalty = corda_options$other_penalty,
    negative_penalty = corda_options$negative_penalty
  )
  direction_model <- parent
  direction_model$lb <- parent$fastcc_original_lb %||% parent$lb
  direction_model$ub <- parent$fastcc_original_ub %||% parent$ub
  target_directions <- rc_prepare_directional_targets(
    direction_model, core, target_direction = target_direction
  )
  parent_diagnostics <- .rc_directional_feasibility(
    parent, target_directions,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  parent_feasible_targets <- parent_diagnostics[
    parent_diagnostics$feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  initial <- .rc_subset_gem(parent, classes$biological)
  initial_diagnostics <- .rc_directional_feasibility(
    initial, parent_feasible_targets,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  names(initial_diagnostics)[names(initial_diagnostics) == "feasible"] <-
    "initial_feasible"
  names(initial_diagnostics)[names(initial_diagnostics) == "vmax"] <-
    "initial_vmax"
  names(initial_diagnostics)[names(initial_diagnostics) == "solver_status"] <-
    "initial_solver_status"
  blocked <- initial_diagnostics[
    !initial_diagnostics$initial_feasible,
    c("reaction_id", "target_direction"),
    drop = FALSE
  ]
  selected_support <- character()
  completion_iterations <- list()
  for (direction in c("forward", "reverse")) {
    task <- blocked[blocked$target_direction == direction, , drop = FALSE]
    if (!nrow(task)) next
    completed <- .rc_corda_complete_direction(
      parent = parent,
      biological_reactions = classes$biological,
      selected_support = selected_support,
      targets = task,
      direction = direction,
      epsilon = fastcore_epsilon,
      solver = solver,
      time_limit = time_limit,
      max_support_reactions = max_support_reactions,
      support_costs = support_costs
    )
    selected_support <- completed$support
    if (nrow(completed$iterations)) {
      completion_iterations[[direction]] <- completed$iterations
    }
  }
  final_reactions <- union(classes$biological, selected_support)
  final <- .rc_subset_gem(parent, final_reactions)
  final_diagnostics <- .rc_directional_feasibility(
    final, parent_feasible_targets,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = fastcore_epsilon
  )
  names(final_diagnostics)[names(final_diagnostics) == "feasible"] <-
    "final_feasible"
  names(final_diagnostics)[names(final_diagnostics) == "vmax"] <-
    "final_vmax"
  names(final_diagnostics)[names(final_diagnostics) == "solver_status"] <-
    "final_solver_status"
  diagnostics <- merge(
    parent_diagnostics, initial_diagnostics,
    by = c("reaction_id", "target_direction"), all.x = TRUE, sort = FALSE
  )
  diagnostics <- merge(
    diagnostics, final_diagnostics,
    by = c("reaction_id", "target_direction"), all.x = TRUE, sort = FALSE
  )
  diagnostics$completion_status <- ifelse(
    diagnostics$target_direction == "none",
    "no_allowed_direction",
    ifelse(
      !diagnostics$feasible,
      "parent_blocked",
      ifelse(
        diagnostics$initial_feasible %in% TRUE,
        "already_feasible",
        ifelse(
          diagnostics$final_feasible %in% TRUE,
          "corda_like_weighted_completed",
          "unresolved"
        )
      )
    )
  )
  failed <- diagnostics$feasible %in% TRUE &
    !(diagnostics$final_feasible %in% TRUE)
  if (isTRUE(strict) && any(failed)) {
    bad <- paste(
      paste(
        diagnostics$reaction_id[failed],
        diagnostics$target_direction[failed],
        sep = ":"
      ),
      collapse = ", "
    )
    stop(
      "Cell-type CORDA-like union-GEM completion failed for parent-feasible ",
      "targets in `", cell_type, "`: ", bad,
      call. = FALSE
    )
  }
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
  meta$celltype_fastcore_support <- reaction_id %in% selected_support
  meta$celltype_corda_support <- reaction_id %in% selected_support
  meta$support_only <- reaction_id %in% selected_support &
    !reaction_id %in% classes$biological
  meta$corda_evidence_class <- unname(classes$evidence_class[reaction_id])
  meta$corda_evidence_score <- unname(classes$evidence_score[reaction_id])
  meta$corda_rna_percentile <- suppressWarnings(as.numeric(
    evidence$rna_percentile
  ))
  meta$corda_multiome_percentile <- suppressWarnings(as.numeric(
    evidence$multiome_percentile
  ))
  meta$corda_regulatory_support <- suppressWarnings(as.numeric(
    evidence$regulatory_support
  ))
  meta$corda_support_penalty <- unname(support_costs[reaction_id])
  final$reaction_meta <- meta
  final$sample_id <- cell_type
  final$grn_module_id <- paste0("CELLTYPE_MEDIUM_UNION_GEM::", cell_type)
  final$cell_type <- cell_type
  final$target_directions <- parent_feasible_targets
  final$closure_diagnostics <- diagnostics
  final$completion_iterations <- if (length(completion_iterations)) {
    do.call(rbind, completion_iterations)
  } else {
    data.frame()
  }
  final$corda_reaction_evidence <- reaction_evidence
  n_no_direction <- sum(diagnostics$completion_status == "no_allowed_direction")
  n_parent_blocked <- sum(diagnostics$completion_status == "parent_blocked")
  final$target_status <- if (any(failed)) {
    "structurally_infeasible"
  } else if (nrow(diagnostics) > 0L && n_no_direction == nrow(diagnostics)) {
    "no_allowed_direction"
  } else if (n_no_direction > 0L) {
    "partial_no_allowed_direction"
  } else if (nrow(diagnostics) > 0L &&
             n_parent_blocked == nrow(diagnostics)) {
    "parent_blocked"
  } else if (n_parent_blocked > 0L) {
    "partial_parent_blocked"
  } else {
    "ok"
  }
  final$is_union_gem <- TRUE
  final$union_gem_scope <-
    "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
  final$build_params <- list(
    strategy = "celltype_medium_corda_like_evidence_max",
    cell_type = cell_type,
    algorithm =
      "retain_all_HC_and_MC_then_weighted_add_only_flux_consistent_completion",
    completion_stage =
      "celltype_specific_corda_like_after_condition_module_union",
    evidence_schema = "regcompass_corda_like_reaction_evidence_v1",
    n_celltype_biological_reactions = length(classes$biological),
    n_celltype_fastcore_support_reactions = length(selected_support),
    n_high_confidence_reactions = length(classes$hc),
    n_module_medium_confidence_reactions = length(classes$mc_module),
    n_evidence_medium_confidence_reactions = length(classes$mc_evidence),
    n_other_reactions = length(classes$ot),
    n_negative_confidence_reactions = length(classes$nc),
    n_corda_support_reactions = length(selected_support),
    n_other_support_reactions = sum(selected_support %in% classes$ot),
    n_negative_support_reactions = sum(selected_support %in% classes$nc),
    n_fastcc_consistent_parent_reactions = length(
      parent$fastcc_consistent_reactions %||% colnames(parent$S)
    ),
    n_fastcc_inconsistent_parent_reactions = length(
      parent$fastcc_inconsistent_reactions %||% character()
    ),
    fastcore_epsilon = fastcore_epsilon,
    target_direction = target_direction,
    high_confidence_reactions = classes$hc,
    medium_confidence_reactions = classes$medium_confidence,
    selected_support_reactions = selected_support,
    forbidden_roles = c("demand", "sink", "artificial_support"),
    strict = strict,
    corda_options = corda_options
  )
  final
}
