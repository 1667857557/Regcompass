.rc_bind_frames_fill <- function(x) {
  x <- x[vapply(
    x,
    function(z) is.data.frame(z) && nrow(z) > 0L,
    logical(1)
  )]
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(columns, colnames(z))
    for (column in missing) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.rc_phase_bpparam <- function(
    workers = NULL,
    backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  if (identical(backend, "serial")) return(FALSE)
  param <- rc_default_bpparam(workers = workers, backend = backend)
  param %||% FALSE
}

.rc_release_bpparam <- function(param) {
  if (!identical(param, FALSE) && !is.null(param) &&
      requireNamespace("BiocParallel", quietly = TRUE)) {
    try(BiocParallel::bpstop(param), silent = TRUE)
  }
  invisible(gc(verbose = FALSE))
}

.rc_metacell_logcpm <- function(
    counts, scale_factor = 1e6, library_size = NULL) {
  counts <- .rc_as_dgCMatrix(counts)
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L ||
      !is.finite(scale_factor) || scale_factor <= 0) {
    stop("`scale_factor` must be one positive finite number.", call. = FALSE)
  }
  normalization_scope <- "input_matrix_library_size"
  if (is.null(library_size)) {
    library_size <- Matrix::colSums(counts)
  } else {
    normalization_scope <- "full_transcriptome_library_size_before_gpr_filter"
    if (!is.null(names(library_size))) {
      missing <- setdiff(colnames(counts), names(library_size))
      if (length(missing)) {
        stop(
          "`library_size` is missing metacells: ",
          paste(utils::head(missing, 10L), collapse = ", "),
          call. = FALSE
        )
      }
      library_size <- library_size[colnames(counts)]
    }
  }
  library_size <- as.numeric(library_size)
  if (length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop(
      "`library_size` must contain one positive finite value per metacell.",
      call. = FALSE
    )
  }
  scaled <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  answer <- log1p(scaled)
  dimnames(answer) <- dimnames(counts)
  attr(answer, "normalization_scope") <- normalization_scope
  attr(answer, "library_size") <- stats::setNames(
    library_size,
    colnames(counts)
  )
  answer
}

.rc_local_fastcore_rows <- function(model, membership, core_ids) {
  reactions <- colnames(model$S)
  template <- membership[rep(1L, length(reactions)), , drop = FALSE]
  template$reaction_id <- reactions
  matched <- match(reactions, as.character(membership$reaction_id))
  existing <- !is.na(matched)
  if (any(existing)) {
    template[existing, colnames(membership)] <- membership[
      matched[existing],
      colnames(membership),
      drop = FALSE
    ]
  }
  support <- character()
  if (!is.null(model$reaction_meta) &&
      all(c("reaction_id", "fastcore_support") %in%
          colnames(model$reaction_meta))) {
    support <- as.character(model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support %in% TRUE
    ])
  }
  template$is_core <- reactions %in% core_ids
  template$biological_meta_module_member <- reactions %in%
    as.character(membership$reaction_id)
  template$local_fastcore_support <- reactions %in% support
  previous_stage <- if ("inclusion_stage" %in% colnames(template)) {
    as.character(template$inclusion_stage)
  } else {
    rep(NA_character_, nrow(template))
  }
  previous_stage[is.na(previous_stage) | !nzchar(previous_stage)] <-
    "biological_meta_module_member"
  template$inclusion_stage <- ifelse(
    template$local_fastcore_support,
    "local_fastcore_support",
    previous_stage
  )
  template
}

.rc_complete_stratum_meta_modules <- function(
    meta_modules, gem, outdir, local_fastcore_args = list()) {
  defaults <- list(
    enabled = TRUE,
    target_direction = "both",
    solver = "highs",
    time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE,
    save_models = TRUE
  )
  defaults[names(local_fastcore_args)] <- NULL
  args <- c(defaults, local_fastcore_args)
  if (!isTRUE(args$enabled)) {
    return(list(
      completed_reaction_membership = meta_modules$reaction_membership,
      summary = data.frame(),
      diagnostics = data.frame(),
      completion_iterations = data.frame(),
      parent_scope = "disabled"
    ))
  }
  membership <- meta_modules$reaction_membership
  core <- meta_modules$core_gene_reaction
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(membership) ||
      !all(required %in% colnames(membership))) {
    stop(
      "Meta-module reaction membership is incomplete before local FASTCORE.",
      call. = FALSE
    )
  }
  if (!is.data.frame(core) || !all(required %in% colnames(core))) {
    stop(
      "Meta-module core reactions are incomplete before local FASTCORE.",
      call. = FALSE
    )
  }
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  parent <- .rc_fastcore_parent(
    gem,
    medium_table = NULL,
    condition = NULL,
    solver = args$solver,
    time_limit = args$time_limit,
    fastcore_epsilon = args$fastcore_epsilon
  )
  module_keys <- unique(
    membership[, c("sample_id", "module_id"), drop = FALSE]
  )
  completed_rows <- vector("list", nrow(module_keys))
  summaries <- vector("list", nrow(module_keys))
  diagnostics <- vector("list", nrow(module_keys))
  iterations <- vector("list", nrow(module_keys))
  model_dir <- file.path(outdir, "local_fastcore_models")
  if (isTRUE(args$save_models)) {
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  }
  for (index in seq_len(nrow(module_keys))) {
    sample_id <- as.character(module_keys$sample_id[[index]])
    module_id <- as.character(module_keys$module_id[[index]])
    in_module <- as.character(membership$sample_id) == sample_id &
      as.character(membership$module_id) == module_id
    membership_i <- membership[in_module, , drop = FALSE]
    in_core <- as.character(core$sample_id) == sample_id &
      as.character(core$module_id) == module_id
    core_i <- core[in_core, , drop = FALSE]
    if (!nrow(core_i)) {
      stop(
        "No core reactions remain for local meta-module `", module_id, "`.",
        call. = FALSE
      )
    }
    model <- .rc_complete_meta_module(
      gem = gem,
      reaction_membership = membership_i,
      core_reactions = core_i,
      sample_id = sample_id,
      module_id = module_id,
      medium_table = NULL,
      condition = NULL,
      parent_gem = parent,
      target_direction = args$target_direction,
      solver = args$solver,
      time_limit = args$time_limit,
      fastcore_epsilon = args$fastcore_epsilon,
      max_support_reactions = args$max_support_reactions,
      strict = args$strict
    )
    completed_rows[[index]] <- .rc_local_fastcore_rows(
      model,
      membership_i,
      unique(as.character(core_i$reaction_id))
    )
    support_ids <- model$reaction_meta$reaction_id[
      model$reaction_meta$fastcore_support %in% TRUE
    ]
    summaries[[index]] <- data.frame(
      sample_id = sample_id,
      module_id = module_id,
      n_biological_reactions = nrow(membership_i),
      n_core_reactions = length(unique(as.character(core_i$reaction_id))),
      n_local_fastcore_support = length(unique(as.character(support_ids))),
      n_completed_reactions = ncol(model$S),
      target_status = model$target_status,
      parent_scope = "unconstrained_shared_parent_with_fastcc",
      stringsAsFactors = FALSE
    )
    if (nrow(model$closure_diagnostics)) {
      diagnostic_i <- model$closure_diagnostics
      diagnostic_i$sample_id <- sample_id
      diagnostic_i$module_id <- module_id
      diagnostics[[index]] <- diagnostic_i
    }
    if (nrow(model$completion_iterations)) {
      iteration_i <- model$completion_iterations
      iteration_i$sample_id <- sample_id
      iteration_i$module_id <- module_id
      iterations[[index]] <- iteration_i
    }
    if (isTRUE(args$save_models)) {
      safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", module_id)
      saveRDS(model, file.path(model_dir, paste0(safe_id, ".rds")))
    }
  }
  list(
    completed_reaction_membership = .rc_bind_frames_fill(completed_rows),
    summary = .rc_bind_frames_fill(summaries),
    diagnostics = .rc_bind_frames_fill(diagnostics),
    completion_iterations = .rc_bind_frames_fill(iterations),
    parent_scope = "unconstrained_shared_parent_with_fastcc"
  )
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  names_to_merge <- c(
    "sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  out <- lapply(names_to_merge, function(name) {
    .rc_bind_frames_fill(lapply(
      artifacts,
      function(artifact) artifact$grn_meta_modules[[name]]
    ))
  })
  names(out) <- names_to_merge

  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]

  completed <- .rc_bind_frames_fill(lapply(artifacts, function(artifact) {
    artifact$grn_meta_modules$local_completed_reaction_membership %||%
      artifact$grn_meta_modules$reaction_membership
  }))
  if (!length(core_ids) || !nrow(completed)) {
    stop("No completed global meta-module reactions were produced.",
         call. = FALSE)
  }

  biological <- out$reaction_membership
  biological_ids <- unique(as.character(biological$reaction_id))
  completed_ids <- unique(as.character(completed$reaction_id))
  completed_ids <- completed_ids[!is.na(completed_ids) & nzchar(completed_ids)]

  out$biological_reaction_membership <- biological
  out$local_completed_reaction_membership <- completed
  out$local_fastcore_summary <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$grn_meta_modules$local_fastcore_summary
  ))
  out$local_fastcore_diagnostics <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$grn_meta_modules$local_fastcore_diagnostics
  ))
  out$local_fastcore_completion_iterations <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) {
      artifact$grn_meta_modules$local_fastcore_completion_iterations
    }
  ))
  out$global_core_reactions <- data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = core_ids,
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  out$global_reaction_membership <- data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    reaction_id = completed_ids,
    is_core = completed_ids %in% core_ids,
    inclusion_stage = ifelse(
      completed_ids %in% core_ids,
      "global_union_core",
      ifelse(
        completed_ids %in% biological_ids,
        "global_union_biological_member",
        "global_union_local_fastcore_support"
      )
    ),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_global_meta_module_v3"
  out$source_group_ids <- unique(vapply(
    artifacts,
    function(artifact) as.character(artifact$group_id),
    character(1)
  ))
  out$global_union_source <-
    "deduplicated_local_fastcore_completed_meta_modules"
  out
}
