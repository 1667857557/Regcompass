# Common model integration for corrected Python CORDA2 reconstruction.

.rc_corda_core_closure <- function(
    parent, final, core, target_direction, solver, time_limit,
    flux_threshold) {
  direction_model <- parent
  direction_model$lb <- parent$lb
  direction_model$ub <- parent$ub
  requested <- rc_prepare_directional_targets(
    direction_model, core, target_direction = target_direction
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
      "corda2_reaction_retained",
      "corda2_reaction_unresolved"
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

.rc_complete_celltype_medium_corda_gem <- function(
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
  split <- .rc_corda_split_model(
    parent, tolerance = corda_options$flux_tolerance
  )
  reconstruction <- .rc_corda_build_three_stage(
    split = split,
    classes = classes,
    options = corda_options,
    solver = solver,
    time_limit = time_limit
  )
  included <- intersect(reconstruction$included, validated$reactions)
  if (!length(included)) {
    stop("CORDA2 reconstruction retained no reactions.", call. = FALSE)
  }
  final <- .rc_subset_gem(parent, included)
  closure <- .rc_corda_core_closure(
    parent = parent,
    final = final,
    core = core,
    target_direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = corda_options$flux_tolerance
  )
  if (isTRUE(strict) && any(closure$failed)) {
    bad <- paste(
      paste(
        closure$diagnostics$reaction_id[closure$failed],
        closure$diagnostics$target_direction[closure$failed],
        sep = ":"
      ),
      collapse = ", "
    )
    stop(
      "CORDA2 failed to retain parent-feasible core directions in `",
      cell_type, "`: ", bad,
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
  meta$celltype_fastcore_support <- FALSE
  meta$celltype_corda2_included <- TRUE
  meta$celltype_corda_included <- TRUE
  meta$support_only <- !reaction_id %in% module_reactions
  meta$corda2_initial_confidence <- unname(
    classes$initial_confidence[reaction_id]
  )
  meta$corda2_final_confidence <- unname(
    reconstruction$final_confidence[reaction_id]
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
  # Legacy aliases allow draft objects to be inspected with older notebooks.
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
  final$corda_execution <- reconstruction$execution
  final$corda_stage2_nc_support_pairs <-
    reconstruction$stage2_nc_support_pairs
  final$corda_stage2_nc_support_count <-
    reconstruction$stage2_nc_support_count
  final$corda_reconstruction <- reconstruction
  final$target_status <- if (any(closure$failed)) {
    "structurally_infeasible"
  } else if (!nrow(closure$feasible_targets)) {
    "parent_blocked"
  } else {
    "ok"
  }
  final$is_union_gem <- TRUE
  final$union_gem_scope <-
    "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
  initial <- classes$initial_confidence
  final$build_params <- list(
    strategy = "celltype_medium_corrected_python_corda2",
    cell_type = cell_type,
    algorithm = reconstruction$algorithm,
    completion_stage = "corrected_python_CORDA2_after_confidence_mapping",
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
      "Python_CORDA2_directional_confidence_then_restore_original_bounds",
    fastcore_epsilon_used = FALSE,
    max_support_reactions_ignored_for_corda2 = max_support_reactions,
    strict = strict,
    corda_options = corda_options,
    included_reactions = included,
    stage_update_policy = reconstruction$stage_update_policy,
    python_reference_commit = reconstruction$python_reference_commit
  )
  final
}

# Draft-branch aliases; the canonical public option is `corda2`.
.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
