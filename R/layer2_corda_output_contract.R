# Validate original CORDA2 union models and compact diagnostics.

.rc_corda_empty_task_table <- function() {
  data.frame(
    variable_id = character(), reaction_id = character(),
    direction = character(), stage = character(), replicate = integer(),
    kind = character(), status = character(), target_flux = numeric(),
    vmax = numeric(), objective = numeric(), backend = character(),
    solver_message = character(), opposite_direction_blocked = character(),
    n_associated = integer(), corda2_n_solves = integer(),
    associated = character(), task_key = character(),
    stringsAsFactors = FALSE
  )
}

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
    stop("CORDA2 task diagnostics are missing required columns.",
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

.rc_corda_task_keys <- function(tab) {
  required <- c("stage", "kind", "variable_id", "replicate")
  if (!is.data.frame(tab) || !all(required %in% colnames(tab))) {
    stop("CORDA2 task diagnostics cannot construct stable task keys.",
         call. = FALSE)
  }
  if (!nrow(tab)) return(character())
  paste(
    as.character(tab$stage), as.character(tab$kind),
    utils::URLencode(as.character(tab$variable_id), reserved = TRUE),
    as.integer(tab$replicate), sep = "::"
  )
}

.rc_corda2_reaction_from_variable <- function(variable_id) {
  sub("_CORDA_rev_rxn$", "", as.character(variable_id))
}

.rc_corda_normalize_associations <- function(tab) {
  empty <- data.frame(
    task_key = character(), associated_variable_id = character(),
    associated_reaction_id = character(), stringsAsFactors = FALSE
  )
  if (!is.data.frame(tab) || !nrow(tab)) return(empty)
  if (!"associated" %in% colnames(tab)) {
    stop("CORDA2 task diagnostics have no association payload.",
         call. = FALSE)
  }
  task_key <- .rc_corda_task_keys(tab)
  split_value <- strsplit(as.character(tab$associated), ";", fixed = TRUE)
  rows <- lapply(seq_along(split_value), function(i) {
    variable <- unique(split_value[[i]])
    variable <- variable[!is.na(variable) & nzchar(variable)]
    if (!length(variable)) return(NULL)
    variable <- sort(variable)
    data.frame(
      task_key = rep(task_key[[i]], length(variable)),
      associated_variable_id = variable,
      associated_reaction_id = .rc_corda2_reaction_from_variable(variable),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  answer <- unique(do.call(rbind, rows))
  rownames(answer) <- NULL
  answer[order(
    answer$task_key, answer$associated_reaction_id,
    answer$associated_variable_id
  ), , drop = FALSE]
}

.rc_validate_corda_union_model <- function(model, cell_type) {
  validated <- rc_validate_gem(model)
  if (!isTRUE(model$is_union_gem) ||
      !identical(as.character(model$cell_type), as.character(cell_type)) ||
      !identical(
        as.character(model$union_gem_scope),
        "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
      )) {
    stop("CORDA2 union-model provenance is incomplete.", call. = FALSE)
  }

  core <- unique(trimws(as.character(
    model$required_core_reactions %||% character()
  )))
  core <- core[!is.na(core) & nzchar(core)]
  if (!length(core)) {
    stop("CORDA2 union model has no required core-reaction contract.",
         call. = FALSE)
  }
  missing_core <- setdiff(core, validated$reactions)
  if (length(missing_core)) {
    stop(
      "CORDA2 union model dropped required core reactions: ",
      paste(utils::head(missing_core, 10L), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  targets <- as.data.frame(model$target_directions)
  required <- c("reaction_id", "target_direction")
  if (!all(required %in% colnames(targets))) {
    stop("CORDA2 scoring target directions are incomplete.", call. = FALSE)
  }
  if (nrow(targets)) {
    reaction <- as.character(targets$reaction_id)
    direction <- as.character(targets$target_direction)
    if (anyDuplicated(targets[required]) ||
        any(!reaction %in% validated$reactions) ||
        any(!reaction %in% core) ||
        any(!direction %in% c("forward", "reverse"))) {
      stop(
        "CORDA2 scoring targets must be core reactions in the final GEM.",
        call. = FALSE
      )
    }
    index <- match(reaction, validated$reactions)
    direction_allowed <- ifelse(
      direction == "forward", validated$ub[index] > 0,
      validated$lb[index] < 0
    )
    if (any(!direction_allowed)) {
      stop("CORDA2 scoring targets contain a disallowed final direction.",
           call. = FALSE)
    }
  }
  meta <- model$reaction_meta
  if (!is.data.frame(meta) || !"reaction_id" %in% colnames(meta) ||
      !identical(as.character(meta$reaction_id), validated$reactions)) {
    stop("CORDA2 reaction metadata do not align to the final GEM.",
         call. = FALSE)
  }
  build <- model$build_params
  if (!is.list(build) ||
      !identical(
        as.character(build$algorithm),
        "schultzdre_MATLAB_CORDA2_original_semantics"
      ) ||
      !identical(
        as.character(build$stage_update_policy),
        "original_matlab_directional_order"
      ) ||
      !identical(
        as.character(build$core_retention_policy),
        "immutable_structural_backbone"
      ) ||
      !identical(build$post_reconstruction_closure_lp, FALSE) ||
      !identical(
        as.integer(build$n_retained_core_reactions),
        as.integer(length(core))
      )) {
    stop(
      "CORDA2 build parameters do not satisfy the core-retention/scoring contract.",
      call. = FALSE
    )
  }
  task_tab <- model$corda_task_diagnostics
  edge_tab <- model$corda_association_edges
  if (!is.data.frame(task_tab) || !"task_key" %in% colnames(task_tab) ||
      anyDuplicated(task_tab$task_key)) {
    stop("CORDA2 compact task diagnostics have invalid task keys.",
         call. = FALSE)
  }
  edge_required <- c(
    "task_key", "associated_variable_id", "associated_reaction_id"
  )
  if (!is.data.frame(edge_tab) ||
      !all(edge_required %in% colnames(edge_tab)) ||
      any(!edge_tab$task_key %in% task_tab$task_key) ||
      any(.rc_corda2_reaction_from_variable(
        edge_tab$associated_variable_id
      ) != edge_tab$associated_reaction_id)) {
    stop("CORDA2 normalized association edges do not match tasks.",
         call. = FALSE)
  }
  expected <- sum(as.integer(task_tab$n_associated), na.rm = TRUE)
  if (nrow(edge_tab) != expected) {
    stop("CORDA2 normalized association count does not match tasks.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.rc_finalize_corda_union_model <- function(model, cell_type) {
  task_tab <- model$corda_task_diagnostics
  if (!is.data.frame(task_tab) || !nrow(task_tab)) {
    task_tab <- .rc_corda_empty_task_table()
  }
  association_edges <- .rc_corda_normalize_associations(task_tab)
  task_tab$task_key <- .rc_corda_task_keys(task_tab)
  if ("associated" %in% colnames(task_tab)) task_tab$associated <- NULL
  model$corda_task_summary <- .rc_corda_task_summary(task_tab)
  model$corda_task_diagnostics <- task_tab
  model$corda_association_edges <- association_edges

  reconstruction <- model$corda_reconstruction
  if (is.list(reconstruction)) {
    reconstruction$task_diagnostics <- NULL
    model$corda_reconstruction <- reconstruction
  }
  model$build_params$diagnostic_storage <- paste(
    "task metadata stored once in corda_task_diagnostics; directional",
    "associations normalized in corda_association_edges"
  )
  model$build_params$association_edge_schema <- c(
    "task_key", "associated_variable_id", "associated_reaction_id"
  )
  .rc_validate_corda_union_model(model, cell_type = cell_type)
  model
}
