# Original three-stage CORDA reconstruction for one cell type and medium.

.rc_corda_promote <- function(confidence, reactions, stage, inclusion_stage) {
  reactions <- intersect(unique(as.character(reactions)), names(confidence))
  if (length(reactions)) {
    confidence[reactions] <- "RE"
    missing_stage <- is.na(inclusion_stage[reactions]) |
      !nzchar(inclusion_stage[reactions])
    inclusion_stage[reactions[missing_stage]] <- stage
  }
  list(confidence = confidence, inclusion_stage = inclusion_stage)
}

.rc_corda_nc_support_pairs <- function(results, remaining_nc) {
  rows <- lapply(results, function(result) {
    associated <- intersect(result$associated, remaining_nc)
    if (!length(associated)) return(NULL)
    data.frame(
      mc_reaction = as.character(result$task$reaction_id[[1L]]),
      nc_reaction = associated,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      mc_reaction = character(), nc_reaction = character(),
      stringsAsFactors = FALSE
    ))
  }
  unique(do.call(rbind, rows))
}

.rc_corda_build_three_stage <- function(
    split, classes, options, solver, time_limit) {
  confidence <- classes$confidence
  inclusion_stage <- stats::setNames(
    rep(NA_character_, length(confidence)), names(confidence)
  )
  promoted <- .rc_corda_promote(
    confidence, classes$hc, "initial_HC", inclusion_stage
  )
  confidence <- promoted$confidence
  inclusion_stage <- promoted$inclusion_stage
  execution <- list()
  task_tables <- list()

  stage1_tasks <- .rc_corda_make_tasks(
    split, classes$hc,
    stage = "stage1_hc_dependencies",
    n = options$n,
    kind = "dependency"
  )
  stage1 <- .rc_corda_run_tasks(
    split, stage1_tasks, confidence, options,
    solver = solver, time_limit = time_limit
  )
  execution$stage1 <- stage1$execution
  task_tables$stage1 <- .rc_corda_results_table(stage1$results)
  stage1_associated <- .rc_corda_associated(
    stage1$results, allowed = union(classes$mc, classes$nc)
  )
  promoted <- .rc_corda_promote(
    confidence, stage1_associated,
    "stage1_associated_with_HC", inclusion_stage
  )
  confidence <- promoted$confidence
  inclusion_stage <- promoted$inclusion_stage

  remaining_mc <- names(confidence)[grepl("^MC", confidence)]
  remaining_nc <- names(confidence)[confidence == "NC"]
  stage2_tasks <- .rc_corda_make_tasks(
    split, remaining_mc,
    stage = "stage2_mc_nc_support",
    n = options$n,
    kind = "dependency"
  )
  stage2 <- .rc_corda_run_tasks(
    split, stage2_tasks, confidence, options,
    solver = solver, time_limit = time_limit
  )
  execution$stage2_dependency <- stage2$execution
  task_tables$stage2_dependency <- .rc_corda_results_table(stage2$results)
  support_pairs <- .rc_corda_nc_support_pairs(
    stage2$results, remaining_nc
  )
  nc_support_count <- stats::setNames(rep(0L, length(remaining_nc)), remaining_nc)
  if (nrow(support_pairs)) {
    observed <- table(support_pairs$nc_reaction)
    nc_support_count[names(observed)] <- as.integer(observed)
  }
  shared_nc <- names(nc_support_count)[nc_support_count >= options$p]
  promoted <- .rc_corda_promote(
    confidence, shared_nc,
    "stage2_NC_supports_at_least_p_MC", inclusion_stage
  )
  confidence <- promoted$confidence
  inclusion_stage <- promoted$inclusion_stage

  remaining_nc <- names(confidence)[confidence == "NC"]
  split_after_nc_block <- .rc_corda_block_reactions(split, remaining_nc)
  remaining_mc <- names(confidence)[grepl("^MC", confidence)]
  stage2_feasibility_tasks <- .rc_corda_make_tasks(
    split_after_nc_block, remaining_mc,
    stage = "stage2_remaining_mc_feasibility",
    n = 1L,
    kind = "feasibility"
  )
  stage2_feasibility <- .rc_corda_run_tasks(
    split_after_nc_block, stage2_feasibility_tasks,
    confidence, options,
    solver = solver, time_limit = time_limit
  )
  execution$stage2_feasibility <- stage2_feasibility$execution
  task_tables$stage2_feasibility <- .rc_corda_results_table(
    stage2_feasibility$results
  )
  feasible_mc <- unique(vapply(
    stage2_feasibility$results[
      vapply(stage2_feasibility$results, function(result) {
        identical(result$status, "optimal") &&
          is.finite(result$target_flux) &&
          result$target_flux >= options$epsilon
      }, logical(1))
    ],
    function(result) as.character(result$task$reaction_id[[1L]]),
    character(1)
  ))
  promoted <- .rc_corda_promote(
    confidence, feasible_mc,
    "stage2_MC_feasible_after_remaining_NC_block", inclusion_stage
  )
  confidence <- promoted$confidence
  inclusion_stage <- promoted$inclusion_stage

  blocked_stage3 <- names(confidence)[
    confidence == "NC" | grepl("^MC", confidence)
  ]
  split_stage3 <- .rc_corda_block_reactions(split, blocked_stage3)
  re_before_stage3 <- names(confidence)[confidence == "RE"]
  stage3_tasks <- .rc_corda_make_tasks(
    split_stage3, re_before_stage3,
    stage = "stage3_re_ot_dependencies",
    n = options$n,
    kind = "dependency"
  )
  stage3 <- .rc_corda_run_tasks(
    split_stage3, stage3_tasks, confidence, options,
    solver = solver, time_limit = time_limit
  )
  execution$stage3 <- stage3$execution
  task_tables$stage3 <- .rc_corda_results_table(stage3$results)
  stage3_ot <- .rc_corda_associated(
    stage3$results,
    allowed = names(confidence)[confidence == "OT"]
  )
  promoted <- .rc_corda_promote(
    confidence, stage3_ot,
    "stage3_associated_OT", inclusion_stage
  )
  confidence <- promoted$confidence
  inclusion_stage <- promoted$inclusion_stage

  list(
    included = names(confidence)[confidence == "RE"],
    final_confidence = confidence,
    inclusion_stage = inclusion_stage,
    stage1_associated = stage1_associated,
    stage2_nc_support_pairs = support_pairs,
    stage2_nc_support_count = nc_support_count,
    stage2_promoted_nc = shared_nc,
    stage2_promoted_mc = feasible_mc,
    stage3_associated_ot = stage3_ot,
    blocked_after_stage2 = remaining_nc,
    blocked_before_stage3 = blocked_stage3,
    task_diagnostics = .rc_bind_frames_fill(task_tables),
    execution = execution,
    algorithm = "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment",
    stage_update_policy = "barrier_then_union_order_independent"
  )
}

.rc_corda_core_closure <- function(
    parent, final, core, target_direction, solver, time_limit,
    flux_threshold) {
  direction_model <- parent
  direction_model$lb <- parent$fastcc_original_lb %||% parent$lb
  direction_model$ub <- parent$fastcc_original_ub %||% parent$ub
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
      "corda_retained",
      "corda_unresolved"
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
    stop("CORDA requires valid module and core reactions.", call. = FALSE)
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
    stop("CORDA reconstruction retained no reactions.", call. = FALSE)
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
      "CORDA failed to retain parent-feasible HC directions in `",
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
  meta$celltype_corda_included <- TRUE
  meta$support_only <- !reaction_id %in% module_reactions
  meta$corda_initial_confidence <- unname(
    classes$initial_confidence[reaction_id]
  )
  meta$corda_final_confidence <- unname(
    reconstruction$final_confidence[reaction_id]
  )
  meta$corda_inclusion_stage <- unname(
    reconstruction$inclusion_stage[reaction_id]
  )
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
  final$reaction_meta <- meta
  final$sample_id <- cell_type
  final$grn_module_id <- paste0("CELLTYPE_MEDIUM_UNION_GEM::", cell_type)
  final$cell_type <- cell_type
  final$target_directions <- closure$feasible_targets
  final$closure_diagnostics <- closure$diagnostics
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
    strategy = "celltype_medium_original_corda",
    cell_type = cell_type,
    algorithm = reconstruction$algorithm,
    completion_stage = "original_CORDA_three_stage_after_confidence_mapping",
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
    n_fastcc_consistent_parent_reactions = length(
      parent$fastcc_consistent_reactions %||% colnames(parent$S)
    ),
    n_fastcc_inconsistent_parent_reactions = length(
      parent$fastcc_inconsistent_reactions %||% character()
    ),
    scoring_target_direction = target_direction,
    reconstruction_direction_policy = "all_parent_allowed_directions",
    fastcc_epsilon = fastcore_epsilon,
    max_support_reactions_ignored_for_corda = max_support_reactions,
    strict = strict,
    corda_options = corda_options,
    included_reactions = included,
    stage_update_policy = reconstruction$stage_update_policy
  )
  final
}

# Draft-branch compatibility alias; canonical implementation is `corda`.
.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
