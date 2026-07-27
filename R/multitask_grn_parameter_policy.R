# Canonical parameter and observability policy for direct condition-theta GRNs.
# Loaded after the direct fitter and inference functions.

.rc_run_celltype_multitask_grns_pre_policy <- .rc_run_celltype_multitask_grns

.rc_multitask_grn_defaults <- function() {
  list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 1,
    lambda_rule = "lambda.1se",
    nfolds = 5L,
    n_bootstrap = 100L,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_abs_effect = 0,
    min_cv_rsq = 0,
    min_bootstrap_success_fraction = 0.8,
    zero_tolerance = 1e-8,
    candidate_screen_threshold = 0,
    max_edges_per_target = Inf,
    min_detected_cells_per_condition = 10L,
    min_detection_fraction_per_condition = 0.01,
    seed = 12345L
  )
}

.rc_validate_multitask_grn_args <- function(args = list()) {
  if (!is.list(args)) {
    stop("`multitask_args` must be a list.", call. = FALSE)
  }
  defaults <- .rc_multitask_grn_defaults()
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown `multitask_args`: ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  out <- modifyList(defaults, args)

  if (!is.numeric(out$alpha) || length(out$alpha) != 1L ||
      !is.finite(out$alpha) || out$alpha <= 0 || out$alpha >= 1) {
    stop(
      paste(
        "`multitask_args$alpha` must be in (0, 1). A positive lasso",
        "component is required for direct condition-specific sparsity, and a",
        "positive ridge component stabilizes correlated TF-peak predictors."
      ),
      call. = FALSE
    )
  }

  for (name in c("global_penalty_factor", "deviation_penalty_factor")) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("`multitask_args$", name, "` must be positive.", call. = FALSE)
    }
  }
  if (!isTRUE(all.equal(
    as.numeric(out$global_penalty_factor),
    as.numeric(out$deviation_penalty_factor),
    tolerance = sqrt(.Machine$double.eps)
  ))) {
    stop(
      paste(
        "`global_penalty_factor` and `deviation_penalty_factor` must be equal.",
        "The direct condition-theta model applies one common elastic-net",
        "penalty to theta[e,c]; the two names remain compatibility aliases."
      ),
      call. = FALSE
    )
  }
  out$global_penalty_factor <- as.numeric(out$global_penalty_factor)
  out$deviation_penalty_factor <- out$global_penalty_factor

  for (name in c("min_selection_frequency", "min_sign_stability")) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1) {
      stop("`multitask_args$", name, "` must be in [0, 1].",
           call. = FALSE)
    }
  }
  value <- out$min_bootstrap_success_fraction
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0 || value > 1) {
    stop(
      "`multitask_args$min_bootstrap_success_fraction` must be in (0, 1].",
      call. = FALSE
    )
  }
  out$min_bootstrap_success_fraction <- as.numeric(value)

  for (name in c("min_abs_effect", "min_cv_rsq", "zero_tolerance")) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      stop("`multitask_args$", name, "` must be non-negative.",
           call. = FALSE)
    }
  }
  if (!is.numeric(out$candidate_screen_threshold) ||
      length(out$candidate_screen_threshold) != 1L ||
      !is.finite(out$candidate_screen_threshold) ||
      !identical(as.numeric(out$candidate_screen_threshold), 0)) {
    stop(
      paste(
        "`multitask_args$candidate_screen_threshold` must be 0. Outcome-based",
        "full-data correlation screening would leak target information into",
        "cross-validation; use the structural Pando design and condition-aware",
        "observability filter instead."
      ),
      call. = FALSE
    )
  }
  out$candidate_screen_threshold <- 0

  max_edges <- out$max_edges_per_target
  if (!is.numeric(max_edges) || length(max_edges) != 1L ||
      is.na(max_edges) || !is.infinite(max_edges) || max_edges < 0) {
    stop(
      paste(
        "`multitask_args$max_edges_per_target` must remain `Inf` in the",
        "canonical model. Pando candidate order is deterministic but is not",
        "an evidence ranking, so a finite top-K would be arbitrary."
      ),
      call. = FALSE
    )
  }
  out$max_edges_per_target <- Inf

  for (name in c(
    "nfolds", "n_bootstrap", "seed", "min_detected_cells_per_condition"
  )) {
    value <- out[[name]]
    minimum <- if (identical(name, "nfolds")) 3L else 1L
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum ||
        abs(value - round(value)) > sqrt(.Machine$double.eps)) {
      stop(
        "`multitask_args$", name, "` must be an integer of at least ",
        minimum, ".", call. = FALSE
      )
    }
    out[[name]] <- as.integer(value)
  }
  fraction <- out$min_detection_fraction_per_condition
  if (!is.numeric(fraction) || length(fraction) != 1L ||
      !is.finite(fraction) || fraction < 0 || fraction > 1) {
    stop(
      "`multitask_args$min_detection_fraction_per_condition` must be in [0, 1].",
      call. = FALSE
    )
  }
  out$min_detection_fraction_per_condition <- as.numeric(fraction)
  out$lambda_rule <- match.arg(
    as.character(out$lambda_rule), c("lambda.1se", "lambda.min")
  )
  out
}

.rc_cv_predictive_gate <- function(cv_rsq, min_cv_rsq = 0) {
  cv_rsq <- suppressWarnings(as.numeric(cv_rsq))
  threshold <- suppressWarnings(as.numeric(min_cv_rsq))
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    stop("`min_cv_rsq` must be one finite non-negative number.", call. = FALSE)
  }
  keep <- is.finite(cv_rsq) & cv_rsq > 0
  if (threshold > 0) keep <- keep & cv_rsq >= threshold
  keep
}

.rc_binomial_wilson_interval <- function(success, total, level = 0.95) {
  success <- suppressWarnings(as.numeric(success))
  total <- suppressWarnings(as.numeric(total))
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- success / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * total)) / total) /
    denominator
  lower <- pmax(0, center - half)
  upper <- pmin(1, center + half)
  invalid <- !is.finite(success) | !is.finite(total) | total <= 0 |
    success < 0 | success > total
  lower[invalid] <- NA_real_
  upper[invalid] <- NA_real_
  list(lower = lower, upper = upper)
}

.rc_fit_multitask_target <- function(
    edges, target, rna, atac, meta, condition_col, args) {
  answer <- .rc_fit_multitask_target_direct(
    edges = edges,
    target = target,
    rna = rna,
    atac = atac,
    meta = meta,
    condition_col = condition_col,
    args = args
  )
  condition <- answer$condition
  if (!is.data.frame(condition) || !nrow(condition)) return(answer)

  cv_ok <- .rc_cv_predictive_gate(condition$cv_rsq, args$min_cv_rsq)
  condition$cv_predictive_above_null <- cv_ok
  condition$active_edge <- condition$active_edge %in% TRUE & cv_ok
  bootstrap_total <- suppressWarnings(as.numeric(
    condition$n_bootstrap_success
  ))
  bootstrap_selected <- condition$selection_frequency * bootstrap_total
  interval <- .rc_binomial_wilson_interval(
    bootstrap_selected, bootstrap_total
  )
  condition$selection_frequency_mc_se <- sqrt(
    condition$selection_frequency * (1 - condition$selection_frequency) /
      bootstrap_total
  )
  condition$selection_frequency_lower_95 <- interval$lower
  condition$selection_frequency_upper_95 <- interval$upper
  condition$sign_agreement_fraction <- ifelse(
    condition$selection_frequency > 0,
    (1 + condition$sign_stability) / 2,
    NA_real_
  )
  answer$condition <- condition

  global <- answer$global
  if (is.data.frame(global) && nrow(global)) {
    global$cv_predictive_above_null <- .rc_cv_predictive_gate(
      global$cv_rsq, args$min_cv_rsq
    )
    answer$global <- global
  }
  diagnostics <- answer$diagnostics
  if (is.data.frame(diagnostics) && nrow(diagnostics)) {
    diagnostics$cv_predictive_above_null <- .rc_cv_predictive_gate(
      diagnostics$cv_rsq, args$min_cv_rsq
    )
    diagnostics$n_active_condition_edges <- sum(
      condition$active_edge %in% TRUE
    )
    diagnostics$n_active_conditions <- length(unique(
      condition$condition[condition$active_edge %in% TRUE]
    ))
    diagnostics$coefficient_parameterization <- "direct_condition_theta"
    diagnostics$theta_penalty_factor <-
      as.numeric(args$global_penalty_factor)
    diagnostics$cv_activation_rule <- if (args$min_cv_rsq > 0) {
      paste0("cv_rsq >= ", args$min_cv_rsq, " and cv_rsq > 0")
    } else {
      "strictly positive out-of-fold cv_rsq"
    }
    answer$diagnostics <- diagnostics
  }
  answer
}

.rc_grn_policy_md5 <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 2, compress = FALSE)
  paste0("md5:", unname(as.character(tools::md5sum(path))))
}

.rc_filter_shared_candidate_observability <- function(
    candidates, rna, atac, meta, condition_col, args,
    design_fingerprint) {
  if (!is.data.frame(candidates) || !nrow(candidates)) {
    stop("Pando produced no structural TF-peak-target candidates.",
         call. = FALSE)
  }
  cells <- colnames(rna)
  if (!identical(cells, colnames(atac))) {
    stop("RNA and ATAC cells must align before candidate filtering.",
         call. = FALSE)
  }
  meta_index <- match(cells, rownames(meta))
  if (anyNA(meta_index)) {
    stop("Candidate-filter metadata do not align with assay cells.",
         call. = FALSE)
  }
  condition <- trimws(as.character(meta[[condition_col]][meta_index]))
  conditions <- sort(unique(condition))
  tf_ids <- unique(as.character(candidates$tf_feature_id))
  peak_ids <- unique(as.character(candidates$atac_feature_id))
  target_ids <- unique(as.character(candidates$target_feature_id))
  missing_rna <- setdiff(unique(c(tf_ids, target_ids)), rownames(rna))
  missing_atac <- setdiff(peak_ids, rownames(atac))
  if (length(missing_rna) || length(missing_atac)) {
    stop("Candidate observability references absent RNA or ATAC features.",
         call. = FALSE)
  }

  observable <- matrix(
    FALSE, nrow = nrow(candidates), ncol = length(conditions),
    dimnames = list(NULL, conditions)
  )
  required_cells <- integer(length(conditions))
  names(required_cells) <- conditions
  for (condition_value in conditions) {
    index <- which(condition == condition_value)
    required <- min(
      length(index),
      max(
        args$min_detected_cells_per_condition,
        ceiling(args$min_detection_fraction_per_condition * length(index))
      )
    )
    required_cells[[condition_value]] <- as.integer(required)
    target_count <- Matrix::rowSums(
      rna[target_ids, index, drop = FALSE] > 0
    )
    target_ok <- as.numeric(target_count[
      match(candidates$target_feature_id, target_ids)
    ]) >= required
    predictor_count <- integer(nrow(candidates))
    chunks <- split(
      seq_len(nrow(candidates)),
      ceiling(seq_len(nrow(candidates)) / 500L)
    )
    for (rows in chunks) {
      tf_present <- as.matrix(
        rna[candidates$tf_feature_id[rows], index, drop = FALSE] > 0
      )
      peak_present <- as.matrix(
        atac[candidates$atac_feature_id[rows], index, drop = FALSE] > 0
      )
      predictor_count[rows] <- rowSums(tf_present & peak_present)
    }
    observable[, condition_value] <-
      predictor_count >= required & target_ok
  }
  keep <- rowSums(observable) > 0L
  if (!any(keep)) {
    stop(
      paste(
        "No structural candidate has an observable TF-RNA x peak-ATAC",
        "predictor and target RNA within any condition. Lower the condition-aware",
        "detection floor only after checking assay depth and cell counts."
      ),
      call. = FALSE
    )
  }
  out <- candidates[keep, , drop = FALSE]
  observable <- observable[keep, , drop = FALSE]
  out$n_observable_conditions <- rowSums(observable)
  out$observable_conditions <- apply(observable, 1L, function(value) {
    paste(colnames(observable)[value], collapse = ";")
  })
  policy <- list(
    min_detected_cells_per_condition =
      args$min_detected_cells_per_condition,
    min_detection_fraction_per_condition =
      args$min_detection_fraction_per_condition,
    required_cells_by_condition = required_cells,
    rule = paste(
      "retain an edge when the TF-RNA x peak-ATAC predictor and target RNA each",
      "meet the same condition-specific detection threshold in at least one",
      "condition"
    )
  )
  out$model_edge_universe_id <- .rc_grn_policy_md5(list(
    pando_design_fingerprint = design_fingerprint,
    edge_ids = sort(as.character(out$edge_id)),
    observability_policy = policy
  ))
  list(
    candidates = out,
    policy = policy,
    n_structural_candidates = nrow(candidates),
    n_observable_candidates = nrow(out),
    n_filtered_candidates = nrow(candidates) - nrow(out),
    model_edge_universe_id = unique(out$model_edge_universe_id)
  )
}

.rc_fit_multitask_celltype_grn <- function(
    design, rna, atac, meta, condition_col,
    multitask_args = list()) {
  args <- .rc_validate_multitask_grn_args(multitask_args)
  Pando::validate_grn_design(design)
  cells <- as.character(design$feature_contract$cell_ids)
  if (!setequal(cells, colnames(rna)) || !setequal(cells, colnames(atac))) {
    stop("Pando design cells and normalized assays do not match.", call. = FALSE)
  }
  rna <- rna[, cells, drop = FALSE]
  atac <- atac[, cells, drop = FALSE]
  meta_index <- match(cells, rownames(meta))
  if (anyNA(meta_index)) {
    stop("Cell metadata do not align to the Pando design.", call. = FALSE)
  }
  meta <- meta[meta_index, , drop = FALSE]
  filtered <- .rc_filter_shared_candidate_observability(
    candidates = design$candidate_edges,
    rna = rna,
    atac = atac,
    meta = meta,
    condition_col = condition_col,
    args = args,
    design_fingerprint = design$design_fingerprint
  )
  candidates <- filtered$candidates
  targets <- unique(as.character(candidates$target))
  fits <- lapply(targets, function(target) {
    .rc_fit_multitask_target(
      edges = candidates[candidates$target == target, , drop = FALSE],
      target = target,
      rna = rna,
      atac = atac,
      meta = meta,
      condition_col = condition_col,
      args = args
    )
  })
  bind <- function(name) {
    .rc_bind_frames_fill(lapply(fits, `[[`, name))
  }
  list(
    global = bind("global"),
    condition = bind("condition"),
    diagnostics = bind("diagnostics"),
    candidates = candidates,
    candidate_observability = filtered,
    params = args
  )
}

.rc_validate_canonical_pando_design_args <- function(args = list()) {
  if (!is.list(args)) {
    stop("`pando_design_args` must be a list.", call. = FALSE)
  }
  defaults <- list(
    peak_to_gene_method = "GREAT",
    upstream = 100000,
    downstream = 0,
    extend = 1000000,
    only_tss = FALSE,
    min_tf_detection = 0,
    min_peak_detection = 0,
    min_target_detection = 0,
    max_edges_per_target = Inf
  )
  out <- modifyList(defaults, args)
  for (name in c(
    "min_tf_detection", "min_peak_detection", "min_target_detection"
  )) {
    if (!is.numeric(out[[name]]) || length(out[[name]]) != 1L ||
        !is.finite(out[[name]]) || out[[name]] < 0 || out[[name]] > 1) {
      stop("`pando_design_args$", name, "` must be in [0, 1].",
           call. = FALSE)
    }
    out[[name]] <- as.numeric(out[[name]])
  }
  max_edges <- out$max_edges_per_target
  if (!is.numeric(max_edges) || length(max_edges) != 1L ||
      is.na(max_edges) || !is.infinite(max_edges) || max_edges < 0) {
    stop(
      paste(
        "`pando_design_args$max_edges_per_target` must remain `Inf` in the",
        "canonical model. Structural candidates are filtered by observability",
        "and regularization, not by an arbitrary deterministic top-K."
      ),
      call. = FALSE
    )
  }
  out$max_edges_per_target <- Inf
  out$peak_to_gene_method <- match.arg(
    as.character(out$peak_to_gene_method), c("GREAT", "Signac")
  )
  for (name in c("upstream", "downstream", "extend")) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      stop("`pando_design_args$", name, "` must be non-negative.",
           call. = FALSE)
    }
  }
  if (!is.logical(out$only_tss) || length(out$only_tss) != 1L ||
      is.na(out$only_tss)) {
    stop("`pando_design_args$only_tss` must be TRUE or FALSE.",
         call. = FALSE)
  }
  out
}

.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 100L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_celltype_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  validated_multitask <- .rc_validate_multitask_grn_args(multitask_args)
  pando_design_args <- .rc_validate_canonical_pando_design_args(
    pando_design_args
  )
  if (is.numeric(min_cells) && length(min_cells) == 1L &&
      is.finite(min_cells) && min_cells < 10L * validated_multitask$nfolds) {
    warning(
      paste0(
        "`min_cells = ", min_cells, "` gives fewer than 10 validation cells ",
        "per condition in each requested CV fold. This is a low-power ",
        "override; the canonical default is 100 cells per condition and cell type."
      ),
      call. = FALSE
    )
  }
  answer <- .rc_run_celltype_multitask_grns_pre_policy(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_cells = min_cells,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_design_args = pando_design_args,
    multitask_args = validated_multitask,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    on_celltype_error = on_celltype_error,
    species = species
  )

  candidates <- answer$tf_peak_gene_candidates
  global <- answer$tf_peak_gene_global
  if (is.data.frame(candidates) && nrow(candidates) &&
      is.data.frame(global) && nrow(global)) {
    key <- paste(candidates[[celltype_col]], candidates$edge_id, sep = "\001")
    observable_columns <- intersect(c(
      celltype_col, "edge_id", "n_observable_conditions",
      "observable_conditions", "model_edge_universe_id"
    ), colnames(global))
    map <- unique(global[, observable_columns, drop = FALSE])
    map_key <- paste(map[[celltype_col]], map$edge_id, sep = "\001")
    index <- match(key, map_key)
    candidates$model_observable <- !is.na(index)
    for (name in setdiff(observable_columns, c(celltype_col, "edge_id"))) {
      candidates[[name]] <- map[[name]][index]
    }
    candidates$n_observable_conditions[
      is.na(candidates$n_observable_conditions)
    ] <- 0L
    answer$tf_peak_gene_candidates <- candidates
  }
  if (is.data.frame(answer$celltype_fit_status) &&
      is.data.frame(answer$tf_peak_gene_candidates)) {
    counts <- stats::aggregate(
      answer$tf_peak_gene_candidates$model_observable %in% TRUE,
      by = list(
        cell_type = answer$tf_peak_gene_candidates[[celltype_col]]
      ), FUN = sum
    )
    names(counts) <- c(celltype_col, "n_observable_candidates")
    index <- match(answer$celltype_fit_status[[celltype_col]],
                   counts[[celltype_col]])
    answer$celltype_fit_status$n_observable_candidates <-
      as.integer(counts$n_observable_candidates[index])
    answer$celltype_fit_status$n_observability_filtered <-
      answer$celltype_fit_status$n_structural_candidates -
      answer$celltype_fit_status$n_observable_candidates
  }
  answer$normalization_policy$pando_design_defaults <- pando_design_args
  answer$normalization_policy$candidate_observability <- list(
    min_detected_cells_per_condition =
      validated_multitask$min_detected_cells_per_condition,
    min_detection_fraction_per_condition =
      validated_multitask$min_detection_fraction_per_condition,
    tf_peak_same_cell_detection = TRUE,
    target_same_condition_detection = TRUE
  )
  answer$normalization_policy$peak_to_gene_policy <- paste(
    "GREAT is the canonical structural peak-to-gene domain rule; it uses",
    "genomic regulatory domains without target-expression correlation",
    "screening before the shared multitask model"
  )
  answer$normalization_policy$coefficient_parameterization <- paste(
    "elastic-net penalties act directly on condition-specific theta[e,c];",
    "beta is the cross-condition mean and delta is derived as theta - beta"
  )
  answer$normalization_policy$penalty_aliases <- paste(
    "global_penalty_factor and deviation_penalty_factor are compatibility",
    "aliases for one common theta penalty and must be equal"
  )
  answer$normalization_policy$bootstrap_precision <- paste(
    "100 full-size condition-stratified bootstrap replicates by default;",
    "selection-frequency Monte Carlo SE and Wilson 95% intervals are reported;",
    "this is bootstrap reproducibility, not formal half-sample PFER control"
  )
  answer$normalization_policy$cv_activation <- paste(
    "active targets require strictly positive out-of-fold R-squared; a positive",
    "user min_cv_rsq adds a stronger floor"
  )
  answer$multitask_args <- validated_multitask
  answer
}
