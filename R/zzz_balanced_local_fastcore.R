.rc_run_regcompass_stratum_v14 <- .rc_run_regcompass_stratum
.rc_merge_stratum_layer1_v14 <- .rc_merge_stratum_layer1
.rc_merge_stratum_meta_modules_v14 <- .rc_merge_stratum_meta_modules

.rc_weighted_quantile <- function(x, weights, probs = 0.5) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[keep]
  weights <- weights[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  order_index <- order(x)
  x <- x[order_index]
  weights <- weights[order_index]
  cumulative <- cumsum(weights) / sum(weights)
  vapply(probs, function(probability) {
    if (!is.finite(probability) || probability < 0 || probability > 1) {
      stop("`probs` values must be finite numbers between 0 and 1.", call. = FALSE)
    }
    if (probability <= 0) return(x[[1L]])
    if (probability >= 1) return(x[[length(x)]])
    index <- which(cumulative >= probability)[[1L]]
    x[[index]]
  }, numeric(1))
}

.rc_equal_sample_weights <- function(sample_ids) {
  sample_ids <- trimws(as.character(sample_ids))
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    stop("Sample IDs used for calibration must be non-missing and non-empty.", call. = FALSE)
  }
  counts <- table(sample_ids)
  weights <- 1 / as.numeric(counts[sample_ids])
  names(weights) <- names(sample_ids)
  weights / sum(weights)
}

.rc_weighted_gene_score <- function(X, weights, min_scale = 0.05, z_clip = 6) {
  X <- as.matrix(X)
  if (length(weights) != ncol(X)) {
    stop("`weights` must contain one value per expression column.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("Gene-score weights must be positive and finite.", call. = FALSE)
  }
  centers <- apply(X, 1L, .rc_weighted_quantile, weights = weights, probs = 0.5)
  mad_sigma <- vapply(seq_len(nrow(X)), function(index) {
    .rc_weighted_quantile(
      abs(X[index, ] - centers[[index]]),
      weights,
      probs = 0.5
    ) * 1.4826
  }, numeric(1))
  iqr_sigma <- vapply(seq_len(nrow(X)), function(index) {
    quantiles <- .rc_weighted_quantile(
      X[index, ],
      weights,
      probs = c(0.25, 0.75)
    )
    diff(quantiles) / 1.349
  }, numeric(1))
  scales <- pmax(mad_sigma, iqr_sigma, min_scale, na.rm = TRUE)
  scales[!is.finite(scales) | scales <= 0] <- min_scale
  z <- sweep(X, 1L, centers, "-")
  z <- sweep(z, 1L, scales, "/")
  z <- pmax(pmin(z, z_clip), -z_clip)
  score <- rc_sigmoid(z)
  dimnames(score) <- dimnames(X)
  score
}

.rc_weighted_q95_calibrate <- function(C_raw, weights, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (length(weights) != ncol(C_raw)) {
    stop("`weights` must contain one value per reaction-capacity column.", call. = FALSE)
  }
  q_values <- apply(
    C_raw,
    1L,
    .rc_weighted_quantile,
    weights = weights,
    probs = 0.95
  )
  C_rel <- sweep(C_raw, 1L, q_values + eps, "/")
  C_rel[C_rel > 1] <- 1
  all_missing <- rowSums(is.finite(C_raw)) == 0L
  if (any(all_missing)) C_rel[all_missing, ] <- NA_real_
  diagnostics <- data.frame(
    reaction_id = rownames(C_raw),
    stratum = "global_sample_balanced",
    n = as.integer(rowSums(is.finite(C_raw))),
    n_global = as.integer(rowSums(is.finite(C_raw))),
    q_stratum = as.numeric(q_values),
    q_stratum_used = as.numeric(q_values),
    q_global = as.numeric(q_values),
    rho_n = 1,
    q_value = as.numeric(q_values),
    quantile_used = 0.95,
    n_finite = as.integer(rowSums(is.finite(C_raw))),
    n_finite_global = as.integer(rowSums(is.finite(C_raw))),
    low_n_flag = rowSums(is.finite(C_raw)) < 20L,
    all_missing_reaction_flag = all_missing,
    sample_balanced = TRUE,
    stringsAsFactors = FALSE
  )
  list(C_rel = C_rel, Q = diagnostics)
}

.rc_normalize_calibration_params <- function(params, sample_col, condition_col,
                                             celltype_col) {
  params <- params %||% list()
  correction <- params$expression_batch_correction %||% "none"
  correction <- match.arg(correction, c("none", "limma"))
  sample_balance <- params$sample_balance %||% TRUE
  if (!is.logical(sample_balance) || length(sample_balance) != 1L || is.na(sample_balance)) {
    stop("`sample_balance` must be TRUE or FALSE.", call. = FALSE)
  }
  technical_batch_cols <- unique(as.character(
    params$technical_batch_cols %||% character()
  ))
  technical_batch_cols <- technical_batch_cols[
    !is.na(technical_batch_cols) & nzchar(technical_batch_cols)
  ]
  preserve_design_cols <- unique(as.character(
    params$preserve_design_cols %||% c(condition_col, celltype_col)
  ))
  preserve_design_cols <- preserve_design_cols[
    !is.na(preserve_design_cols) & nzchar(preserve_design_cols)
  ]
  list(
    sample_balance = sample_balance,
    sample_balance_col = as.character(
      params$sample_balance_col %||% sample_col
    ),
    expression_batch_correction = correction,
    technical_batch_cols = technical_batch_cols,
    preserve_design_cols = preserve_design_cols
  )
}

.rc_apply_limma_batch_correction <- function(X, unit_meta, calibration_params,
                                             sample_col) {
  mode <- calibration_params$expression_batch_correction
  if (identical(mode, "none")) {
    return(list(
      matrix = X,
      diagnostics = data.frame(
        method = "none",
        n_batch_levels = NA_integer_,
        batch_columns = "",
        preserve_columns = paste(
          calibration_params$preserve_design_cols,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    ))
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for expression batch correction.", call. = FALSE)
  }
  batch_cols <- calibration_params$technical_batch_cols
  if (!length(batch_cols)) {
    stop(
      "`technical_batch_cols` is required when `expression_batch_correction = 'limma'`.",
      call. = FALSE
    )
  }
  if (sample_col %in% batch_cols) {
    stop(
      "`sample_id` is a biological replicate and cannot be removed as a technical batch.",
      call. = FALSE
    )
  }
  missing <- setdiff(
    unique(c(batch_cols, calibration_params$preserve_design_cols)),
    colnames(unit_meta)
  )
  if (length(missing)) {
    stop("Batch-correction metadata are missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  batch_frame <- unit_meta[, batch_cols, drop = FALSE]
  if (anyNA(batch_frame)) {
    stop("Technical batch columns contain missing values.", call. = FALSE)
  }
  batch <- do.call(
    interaction,
    c(lapply(batch_frame, as.factor), list(drop = TRUE, lex.order = TRUE))
  )
  n_batch_levels <- nlevels(batch)
  if (n_batch_levels < 2L) {
    warning(
      "Only one technical batch level was detected; limma correction was skipped.",
      call. = FALSE
    )
    return(list(
      matrix = X,
      diagnostics = data.frame(
        method = "limma_skipped_single_level",
        n_batch_levels = n_batch_levels,
        batch_columns = paste(batch_cols, collapse = ";"),
        preserve_columns = paste(
          calibration_params$preserve_design_cols,
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    ))
  }
  preserve_cols <- calibration_params$preserve_design_cols
  preserve_cols <- preserve_cols[vapply(
    unit_meta[, preserve_cols, drop = FALSE],
    function(value) length(unique(value[!is.na(value)])) > 1L,
    logical(1)
  )]
  design <- if (length(preserve_cols)) {
    stats::model.matrix(stats::reformulate(preserve_cols), data = unit_meta)
  } else {
    matrix(1, nrow = nrow(unit_meta), ncol = 1L)
  }
  confounding_data <- unit_meta
  confounding_data$.rc_technical_batch <- batch
  full_terms <- c(preserve_cols, ".rc_technical_batch")
  full_design <- stats::model.matrix(
    stats::reformulate(full_terms),
    data = confounding_data
  )
  if (qr(full_design)$rank < ncol(full_design)) {
    stop(
      "Technical batch is confounded with the preserved biological design; correction is not identifiable.",
      call. = FALSE
    )
  }
  corrected <- limma::removeBatchEffect(
    as.matrix(X),
    batch = batch,
    design = design
  )
  dimnames(corrected) <- dimnames(X)
  list(
    matrix = corrected,
    diagnostics = data.frame(
      method = "limma_removeBatchEffect",
      n_batch_levels = n_batch_levels,
      batch_columns = paste(batch_cols, collapse = ";"),
      preserve_columns = paste(preserve_cols, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
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
      all(c("reaction_id", "fastcore_support") %in% colnames(model$reaction_meta))) {
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

.rc_complete_stratum_meta_modules <- function(meta_modules, gem, outdir,
                                              local_fastcore_args = list()) {
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
  if (!is.data.frame(membership) || !all(required %in% colnames(membership))) {
    stop("Meta-module reaction membership is incomplete before local FASTCORE.", call. = FALSE)
  }
  if (!is.data.frame(core) || !all(required %in% colnames(core))) {
    stop("Meta-module core reactions are incomplete before local FASTCORE.", call. = FALSE)
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
  module_keys <- unique(membership[, c("sample_id", "module_id"), drop = FALSE])
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
      stop("No core reactions remain for local meta-module `", module_id, "`.",
           call. = FALSE)
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

.rc_run_regcompass_stratum <- function(object, group_id, group_cols, gem, outdir,
                                        pfm, genome, fragment_files = NULL,
                                        sample_col = "sample_id",
                                        condition_col = "condition",
                                        celltype_col = "cell_type",
                                        rna_assay = "RNA",
                                        atac_assay = "ATAC",
                                        metacell_args = list(),
                                        layer1_args = list(),
                                        pando_args = list()) {
  result <- .rc_run_regcompass_stratum_v14(
    object = object,
    group_id = group_id,
    group_cols = group_cols,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args
  )
  if (!identical(result$status, "ok")) return(result)
  artifact <- readRDS(result$artifact_file)
  artifact$calibration_params <- .rc_normalize_calibration_params(
    list(
      sample_balance = layer1_args$sample_balance %||% TRUE,
      sample_balance_col = layer1_args$sample_balance_col %||% sample_col,
      expression_batch_correction =
        layer1_args$expression_batch_correction %||% "none",
      technical_batch_cols = layer1_args$technical_batch_cols %||% character(),
      preserve_design_cols = layer1_args$preserve_design_cols %||%
        c(condition_col, celltype_col)
    ),
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  artifact$grn_meta_modules$biological_reaction_membership <-
    artifact$grn_meta_modules$reaction_membership
  local_args <- layer1_args$local_fastcore_args %||% list()
  local_args$enabled <- layer1_args$local_fastcore %||%
    local_args$enabled %||% TRUE
  local_completion <- .rc_complete_stratum_meta_modules(
    artifact$grn_meta_modules,
    gem,
    outdir = file.path(artifact$metacell_dir, "03_local_fastcore"),
    local_fastcore_args = local_args
  )
  artifact$grn_meta_modules$local_completed_reaction_membership <-
    local_completion$completed_reaction_membership
  artifact$grn_meta_modules$local_fastcore_summary <- local_completion$summary
  artifact$grn_meta_modules$local_fastcore_diagnostics <-
    local_completion$diagnostics
  artifact$grn_meta_modules$local_fastcore_completion_iterations <-
    local_completion$completion_iterations
  artifact$grn_meta_modules$local_fastcore_parent_scope <-
    local_completion$parent_scope
  artifact$schema_version <- "regcompass_stratum_v3"
  saveRDS(artifact, result$artifact_file)
  result
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  out <- .rc_merge_stratum_meta_modules_v14(artifacts)
  completed <- .rc_bind_frames_fill(lapply(artifacts, function(artifact) {
    artifact$grn_meta_modules$local_completed_reaction_membership %||%
      artifact$grn_meta_modules$reaction_membership
  }))
  if (!nrow(completed)) {
    stop("No locally completed meta-module reactions were produced.", call. = FALSE)
  }
  out$biological_reaction_membership <- out$reaction_membership
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
    function(artifact) artifact$grn_meta_modules$local_fastcore_completion_iterations
  ))
  core_ids <- unique(as.character(out$global_core_reactions$reaction_id))
  biological_ids <- unique(as.character(out$reaction_membership$reaction_id))
  completed_ids <- unique(as.character(completed$reaction_id))
  completed_ids <- completed_ids[!is.na(completed_ids) & nzchar(completed_ids)]
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
  out$global_union_source <-
    "deduplicated_local_fastcore_completed_meta_modules"
  out
}

.rc_merge_stratum_layer1 <- function(artifacts, gem, single_cell_genes,
                                      sample_col, condition_col, celltype_col) {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  expression_list <- lapply(
    artifacts,
    function(artifact) artifact$layer1$rna_metacell_logcpm
  )
  if (any(vapply(expression_list, is.null, logical(1)))) {
    stop(
      "Every upstream artifact must contain metacell RNA logCPM for global capacity recomputation.",
      call. = FALSE
    )
  }
  expression_list <- lapply(expression_list, function(expression) {
    keep <- tolower(rownames(expression)) %in% gpr_genes
    as.matrix(expression[keep, , drop = FALSE])
  })
  common_genes <- Reduce(intersect, lapply(expression_list, rownames))
  common_genes <- rownames(expression_list[[1L]])[
    rownames(expression_list[[1L]]) %in% common_genes
  ]
  if (!length(common_genes)) {
    stop("No common GPR genes remain across upstream metacell artifacts.", call. = FALSE)
  }
  expression_list <- lapply(
    expression_list,
    function(expression) expression[common_genes, , drop = FALSE]
  )
  rna_logcpm_uncorrected <- .rc_cbind_matrix_union(expression_list)
  unit_meta <- .rc_bind_frames_fill(lapply(
    artifacts,
    function(artifact) artifact$layer1$unit_meta
  ))
  id_col <- if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else {
    NULL
  }
  if (is.null(id_col)) {
    stop("Merged metacell metadata lacks pool_id/metacell_id.", call. = FALSE)
  }
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm_uncorrected), as.character(unit_meta[[id_col]])),
    ,
    drop = FALSE
  ]
  if (anyNA(unit_meta[[id_col]])) {
    stop("Merged metacell metadata are incomplete.", call. = FALSE)
  }
  if (!"pool_id" %in% colnames(unit_meta)) {
    unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  }
  if (!"unit_id" %in% colnames(unit_meta)) unit_meta$unit_id <- unit_meta$pool_id
  unit_meta$stratum_id <- rc_make_stratum_id(
    unit_meta,
    c(condition_col, sample_col, celltype_col)
  )
  parameter_keys <- vapply(artifacts, function(artifact) {
    params <- artifact$capacity_params
    if (is.null(params)) return(NA_character_)
    paste(
      params$promiscuity_mode,
      params$and_method,
      params$tau,
      params$or_method,
      sep = "|"
    )
  }, character(1))
  if (anyNA(parameter_keys) || length(unique(parameter_keys)) != 1L) {
    stop(
      "Upstream artifacts must use one identical Layer 1 capacity parameter set.",
      call. = FALSE
    )
  }
  capacity_params <- artifacts[[1L]]$capacity_params
  calibration_list <- lapply(artifacts, function(artifact) {
    .rc_normalize_calibration_params(
      artifact$calibration_params,
      sample_col = sample_col,
      condition_col = condition_col,
      celltype_col = celltype_col
    )
  })
  calibration_params <- calibration_list[[1L]]
  identical_calibration <- vapply(
    calibration_list,
    identical,
    logical(1),
    y = calibration_params
  )
  if (!all(identical_calibration)) {
    stop("Upstream artifacts use inconsistent global calibration settings.", call. = FALSE)
  }
  correction <- .rc_apply_limma_batch_correction(
    rna_logcpm_uncorrected,
    unit_meta,
    calibration_params,
    sample_col = sample_col
  )
  rna_logcpm <- correction$matrix
  balance_col <- calibration_params$sample_balance_col
  if (!balance_col %in% colnames(unit_meta)) {
    stop("Sample-balance metadata are missing column: ", balance_col, call. = FALSE)
  }
  weights <- if (isTRUE(calibration_params$sample_balance)) {
    .rc_equal_sample_weights(unit_meta[[balance_col]])
  } else {
    rep(1 / ncol(rna_logcpm), ncol(rna_logcpm))
  }
  names(weights) <- colnames(rna_logcpm)
  global_gene_score <- if (isTRUE(calibration_params$sample_balance)) {
    .rc_weighted_gene_score(rna_logcpm, weights)
  } else {
    rc_gene_score(rna_logcpm)
  }
  C_raw <- rc_reaction_capacity(
    parsed,
    global_gene_score,
    promiscuity_mode = capacity_params$promiscuity_mode,
    tau = capacity_params$tau,
    and_method = capacity_params$and_method,
    or_method = capacity_params$or_method,
    BPPARAM = FALSE
  )
  calibrated <- if (isTRUE(calibration_params$sample_balance)) {
    .rc_weighted_q95_calibrate(C_raw, weights)
  } else {
    rc_q95_calibrate(
      C_raw,
      bootstrap = FALSE,
      BPPARAM = FALSE,
      unit_meta = NULL,
      stratum_col = NULL
    )
  }
  confidence_list <- lapply(artifacts, function(artifact) {
    value <- artifact$layer1$reaction_confidence
    if (is.null(value)) {
      stop(
        "Every upstream artifact must contain Pando-derived reaction confidence.",
        call. = FALSE
      )
    }
    as.matrix(value)
  })
  reaction_confidence <- .rc_cbind_matrix_union(confidence_list)
  reaction_confidence <- rc_align_layer2_evidence(
    reaction_confidence,
    rownames(C_raw),
    NA_real_
  )
  missing_confidence_units <- setdiff(
    colnames(C_raw),
    colnames(reaction_confidence)
  )
  if (length(missing_confidence_units)) {
    stop(
      "Reaction confidence is missing metacells after global alignment.",
      call. = FALSE
    )
  }
  reaction_confidence <- reaction_confidence[
    ,
    colnames(C_raw),
    drop = FALSE
  ]
  calibration_scope <- paste(
    if (isTRUE(calibration_params$sample_balance)) {
      "equal_sample_weighted"
    } else {
      "equal_metacell_weighted"
    },
    calibration_params$expression_batch_correction,
    "global_gene_score_and_reaction_q95",
    sep = "_"
  )
  list(
    schema_version = "regcompass_global_layer1_v3",
    C_or_raw = C_raw,
    C_raw = C_raw,
    reaction_capacity_L1 = C_raw,
    C_rel = calibrated$C_rel,
    reaction_confidence = reaction_confidence,
    q95_diagnostics = calibrated$Q,
    capacity_calibration_scope = calibration_scope,
    reaction_confidence_source = "pando_internal_peak_gene_accessibility",
    capacity_params = capacity_params,
    calibration_params = calibration_params,
    sample_balance_weights = weights,
    expression_batch_diagnostics = correction$diagnostics,
    rna_metacell_logcpm_uncorrected = rna_logcpm_uncorrected,
    rna_metacell_logcpm = rna_logcpm,
    global_gene_score = global_gene_score,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, tolower(single_cell_genes)),
    parsed_gpr = parsed,
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "metacell",
    strict_stratum_cols = c(condition_col, sample_col, celltype_col)
  )
}
