# Validate CORDA union models and avoid duplicated large diagnostics.

.rc_complete_celltype_medium_corda_gem_base <-
  .rc_complete_celltype_medium_corda_gem

.rc_corda_task_summary <- function(tab) {
  if (!is.data.frame(tab) || !nrow(tab)) {
    return(data.frame(
      stage = character(), kind = character(), status = character(),
      backend = character(), n_tasks = integer(),
      n_associated_total = integer(), stringsAsFactors = FALSE
    ))
  }
  required <- c("stage", "kind", "status", "backend", "n_associated")
  if (!all(required %in% colnames(tab))) {
    stop("CORDA task diagnostics are missing required columns.",
         call. = FALSE)
  }
  key <- interaction(
    tab$stage, tab$kind, tab$status, tab$backend,
    drop = TRUE, lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(tab)), key), function(index) {
    data.frame(
      stage = as.character(tab$stage[index[[1L]]]),
      kind = as.character(tab$kind[index[[1L]]]),
      status = as.character(tab$status[index[[1L]]]),
      backend = as.character(tab$backend[index[[1L]]]),
      n_tasks = length(index),
      n_associated_total = sum(
        suppressWarnings(as.integer(tab$n_associated[index])), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[do.call(order, answer[c("stage", "kind", "status", "backend")]),
         , drop = FALSE]
}

.rc_validate_corda_union_model <- function(model, cell_type) {
  validated <- rc_validate_gem(model)
  if (!isTRUE(model$is_union_gem) ||
      !identical(as.character(model$cell_type), as.character(cell_type)) ||
      !identical(
        as.character(model$union_gem_scope),
        "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
      )) {
    stop("CORDA union-model provenance is incomplete.", call. = FALSE)
  }
  targets <- as.data.frame(model$target_directions)
  required <- c("reaction_id", "target_direction")
  if (!all(required %in% colnames(targets))) {
    stop("CORDA scoring target directions are incomplete.", call. = FALSE)
  }
  if (nrow(targets)) {
    reaction <- as.character(targets$reaction_id)
    direction <- as.character(targets$target_direction)
    if (anyDuplicated(targets[required]) ||
        any(!reaction %in% validated$reactions) ||
        any(!direction %in% c("forward", "reverse"))) {
      stop("CORDA scoring targets do not match the final GEM.", call. = FALSE)
    }
    index <- match(reaction, validated$reactions)
    direction_allowed <- ifelse(
      direction == "forward",
      validated$ub[index] > 0,
      validated$lb[index] < 0
    )
    if (any(!direction_allowed)) {
      stop("CORDA scoring targets contain a disallowed final direction.",
           call. = FALSE)
    }
  }
  meta <- model$reaction_meta
  if (!is.data.frame(meta) ||
      !"reaction_id" %in% colnames(meta) ||
      !identical(as.character(meta$reaction_id), validated$reactions)) {
    stop("CORDA reaction metadata do not align to the final GEM.",
         call. = FALSE)
  }
  build <- model$build_params
  if (!is.list(build) ||
      !identical(
        as.character(build$algorithm),
        "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment"
      ) ||
      !identical(
        as.character(build$stage_update_policy),
        "barrier_then_union_order_independent"
      )) {
    stop("CORDA build parameters do not identify the published algorithm.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_complete_celltype_medium_corda_gem <- function(...) {
  args <- list(...)
  model <- do.call(.rc_complete_celltype_medium_corda_gem_base, args)
  model$corda_task_summary <- .rc_corda_task_summary(
    model$corda_task_diagnostics
  )
  reconstruction <- model$corda_reconstruction
  if (is.list(reconstruction)) {
    reconstruction$task_diagnostics <- NULL
    reconstruction$execution <- NULL
    reconstruction$stage2_nc_support_pairs <- NULL
    reconstruction$stage2_nc_support_count <- NULL
    model$corda_reconstruction <- reconstruction
  }
  model$build_params$diagnostic_storage <- paste(
    "full task table stored once in corda_task_diagnostics; compact task",
    "summary stored in corda_task_summary; reconstruction stores stage sets"
  )
  .rc_validate_corda_union_model(
    model,
    cell_type = as.character(args$cell_type)
  )
  model
}

.rc_complete_celltype_medium_corda_like_gem <-
  .rc_complete_celltype_medium_corda_gem
