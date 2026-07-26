# Shared-candidate multi-task GRN inference for condition-comparable networks.

.rc_mt_validate_args <- function(args) {
  defaults <- list(
    alpha = 0.5,
    global_penalty_factor = 1,
    deviation_penalty_factor = 2,
    lambda_rule = "lambda.1se",
    nfolds = 5L,
    n_stability = 30L,
    stability_fraction = 0.8,
    min_selection_frequency = 0.7,
    min_sign_stability = 0.8,
    min_abs_effect = 0,
    min_cv_rsq = 0.1,
    interaction_min_scale = 1e-4,
    coefficient_zero_tolerance = 1e-8,
    seed = 12345L
  )
  if (!is.list(args)) stop("`multitask_args` must be a list.", call. = FALSE)
  unknown <- setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop("Unknown `multitask_args`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  answer <- modifyList(defaults, args)
  probability <- function(x, name, open_zero = FALSE) {
    valid <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= if (open_zero) 0 else 0 && x <= 1
    if (!valid || (open_zero && x <= 0)) {
      stop("`", name, "` must be in ", if (open_zero) "(0, 1]" else "[0, 1]",
           ".", call. = FALSE)
    }
  }
  probability(answer$alpha, "alpha")
  if (answer$alpha >= 1) {
    stop(
      "`multitask_args$alpha` must be below 1 so the ridge component makes ",
      "the redundant global/deviation parameterization unique.",
      call. = FALSE
    )
  }
  probability(answer$stability_fraction, "stability_fraction", open_zero = TRUE)
  probability(answer$min_selection_frequency, "min_selection_frequency")
  probability(answer$min_sign_stability, "min_sign_stability")
  positive <- c(
    "global_penalty_factor", "deviation_penalty_factor",
    "interaction_min_scale", "coefficient_zero_tolerance"
  )
  for (name in positive) {
    value <- answer[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0) {
      stop("`multitask_args$", name, "` must be positive.", call. = FALSE)
    }
  }
  if (answer$deviation_penalty_factor <= answer$global_penalty_factor) {
    stop(
      "Condition deviations must be penalized more strongly than the global ",
      "backbone.", call. = FALSE
    )
  }
  nonnegative <- c("min_abs_effect", "min_cv_rsq")
  for (name in nonnegative) {
    value <- answer[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      stop("`multitask_args$", name, "` must be non-negative.",
           call. = FALSE)
    }
  }
  answer$lambda_rule <- match.arg(
    as.character(answer$lambda_rule), c("lambda.1se", "lambda.min")
  )
  integer_fields <- c("nfolds", "n_stability", "seed")
  for (name in integer_fields) {
    value <- answer[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < if (identical(name, "n_stability")) 0 else 1 ||
        abs(value - round(value)) > sqrt(.Machine$double.eps)) {
      stop("`multitask_args$", name, "` must be an integer.", call. = FALSE)
    }
    answer[[name]] <- as.integer(value)
  }
  if (answer$nfolds < 3L) stop("`nfolds` must be at least 3.", call. = FALSE)
  answer
}

.rc_mt_seed <- function(seed, key) {
  values <- utf8ToInt(enc2utf8(as.character(key)))
  offset <- if (length(values)) sum(values * seq_along(values)) else 0
  as.integer((as.double(seed) + offset) %% (.Machine$integer.max - 1L) + 1L)
}

.rc_mt_condition_weights <- function(condition) {
  condition <- factor(condition)
  counts <- table(condition)
  weights <- 1 / (length(counts) * as.numeric(counts[condition]))
  weights / mean(weights)
}

.rc_mt_foldid <- function(condition, nfolds, seed) {
  condition <- factor(condition)
  nfolds <- min(as.integer(nfolds), min(as.integer(table(condition))))
  if (nfolds < 3L) {
    stop("Each condition requires at least three cells for cross-validation.",
         call. = FALSE)
  }
  set.seed(seed)
  foldid <- integer(length(condition))
  for (level in levels(condition)) {
    index <- which(condition == level)
    index <- sample(index, length(index), replace = FALSE)
    foldid[index] <- rep(seq_len(nfolds), length.out = length(index))
  }
  foldid
}

.rc_mt_intercept_matrix <- function(condition, sample = NULL) {
  condition <- as.character(condition)
  if (is.null(sample)) {
    stratum <- condition
    prefix <- "condition"
  } else {
    sample <- as.character(sample)
    stratum <- paste(condition, sample, sep = "::")
    prefix <- "condition_sample"
  }
  levels <- sort(unique(stratum))
  answer <- Matrix::sparseMatrix(
    i = seq_along(stratum),
    j = match(stratum, levels),
    x = 1,
    dims = c(length(stratum), length(levels)),
    dimnames = list(NULL, paste0("I::", prefix, "::", levels))
  )
  list(matrix = answer, stratum = stratum)
}

.rc_mt_balanced_interaction_scale <- function(X, condition, min_scale) {
  condition <- factor(condition)
  variances <- lapply(levels(condition), function(level) {
    selected <- condition == level
    one <- X[selected, , drop = FALSE]
    first <- Matrix::colMeans(one)
    second <- Matrix::colMeans(one ^ 2)
    pmax(as.numeric(second - first ^ 2), 0)
  })
  raw <- sqrt(Reduce("+", variances) / length(variances))
  list(raw = raw, scale = pmax(raw, min_scale))
}

.rc_mt_balanced_tf_reference <- function(tf_matrix, condition) {
  condition <- factor(condition)
  references <- lapply(levels(condition), function(level) {
    as.numeric(Matrix::colMeans(tf_matrix[condition == level, , drop = FALSE]))
  })
  Reduce("+", references) / length(references)
}

.rc_mt_target_design <- function(
    rna, atac, cells, condition, sample, edges, args) {
  tf_matrix <- Matrix::t(rna[edges$tf_feature_id, cells, drop = FALSE])
  peak_matrix <- Matrix::t(atac[edges$atac_feature_id, cells, drop = FALSE])
  X <- tf_matrix * peak_matrix
  colnames(X) <- edges$edge_id
  scaling <- .rc_mt_balanced_interaction_scale(
    X, condition, args$interaction_min_scale
  )
  tf_reference <- .rc_mt_balanced_tf_reference(tf_matrix, condition)
  estimable <- is.finite(scaling$raw) &
    scaling$raw > args$coefficient_zero_tolerance
  edge_metadata <- edges
  edge_metadata$interaction_scale <- scaling$scale
  edge_metadata$raw_interaction_scale <- scaling$raw
  edge_metadata$tf_reference <- tf_reference
  edge_metadata$estimable <- estimable
  if (!any(estimable)) {
    return(list(edge_metadata = edge_metadata, estimable = FALSE))
  }
  edge_metadata <- edge_metadata[estimable, , drop = FALSE]
  X <- X[, estimable, drop = FALSE]
  X <- X %*% Matrix::Diagonal(x = 1 / edge_metadata$interaction_scale)
  colnames(X) <- edge_metadata$edge_id

  intercept <- .rc_mt_intercept_matrix(condition, sample)
  condition_levels <- sort(unique(as.character(condition)))
  blocks <- lapply(condition_levels, function(level) {
    Matrix::Diagonal(x = as.numeric(as.character(condition) == level)) %*% X
  })
  names(blocks) <- condition_levels
  global_names <- paste0("G::", edge_metadata$edge_id)
  colnames(X) <- global_names
  for (level in condition_levels) {
    colnames(blocks[[level]]) <- paste0(
      "D::", level, "::", edge_metadata$edge_id
    )
  }
  design <- do.call(cbind, c(list(intercept$matrix, X), blocks))
  penalty_factor <- c(
    rep(0, ncol(intercept$matrix)),
    rep(args$global_penalty_factor, nrow(edge_metadata)),
    rep(args$deviation_penalty_factor,
        nrow(edge_metadata) * length(condition_levels))
  )
  list(
    design = design,
    penalty_factor = penalty_factor,
    edge_metadata = edge_metadata,
    condition_levels = condition_levels,
    intercept_stratum = intercept$stratum,
    n_intercepts = ncol(intercept$matrix),
    estimable = TRUE
  )
}

.rc_mt_extract_effective <- function(
    coefficients, edge_ids, condition_levels) {
  global <- coefficients[paste0("G::", edge_ids)]
  global[!is.finite(global)] <- 0
  internal_deviation <- vapply(condition_levels, function(condition) {
    value <- coefficients[paste0("D::", condition, "::", edge_ids)]
    value[!is.finite(value)] <- 0
    value
  }, numeric(length(edge_ids)))
  if (is.null(dim(internal_deviation))) {
    internal_deviation <- matrix(
      internal_deviation, nrow = length(edge_ids),
      dimnames = list(edge_ids, condition_levels)
    )
  } else {
    rownames(internal_deviation) <- edge_ids
    colnames(internal_deviation) <- condition_levels
  }
  effective <- sweep(internal_deviation, 1L, global, "+")
  canonical_global <- rowMeans(effective)
  canonical_deviation <- sweep(effective, 1L, canonical_global, "-")
  list(
    global = canonical_global,
    deviation = canonical_deviation,
    effective = effective,
    internal_global = global,
    internal_deviation = internal_deviation
  )
}

.rc_mt_glmnet_coefficients <- function(fit, s, design_names) {
  value <- as.matrix(stats::coef(fit, s = s))
  coefficient <- as.numeric(value[-1L, 1L])
  names(coefficient) <- rownames(value)[-1L]
  answer <- stats::setNames(rep(0, length(design_names)), design_names)
  common <- intersect(names(coefficient), design_names)
  answer[common] <- coefficient[common]
  answer
}

.rc_mt_weighted_rsq <- function(y, prediction, weights) {
  valid <- is.finite(y) & is.finite(prediction) & is.finite(weights) & weights > 0
  if (sum(valid) < 3L) return(NA_real_)
  y <- y[valid]
  prediction <- prediction[valid]
  weights <- weights[valid]
  center <- stats::weighted.mean(y, weights)
  denominator <- sum(weights * (y - center) ^ 2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum(weights * (y - prediction) ^ 2) / denominator
}

.rc_mt_stability_indices <- function(stratum, fraction) {
  rows <- split(seq_along(stratum), stratum)
  sort(unlist(lapply(rows, function(index) {
    size <- max(1L, floor(length(index) * fraction))
    sample(index, min(size, length(index)), replace = FALSE)
  }), use.names = FALSE))
}

.rc_mt_fit_target <- function(
    y, target_design, condition, args, target) {
  if (!isTRUE(target_design$estimable)) {
    stop("No estimable interaction columns for target ", target, ".",
         call. = FALSE)
  }
  design <- target_design$design
  weights <- .rc_mt_condition_weights(condition)
  target_seed <- .rc_mt_seed(args$seed, target)
  foldid <- .rc_mt_foldid(condition, args$nfolds, target_seed)
  cvfit <- glmnet::cv.glmnet(
    x = design,
    y = as.numeric(y),
    family = "gaussian",
    weights = weights,
    foldid = foldid,
    alpha = args$alpha,
    penalty.factor = target_design$penalty_factor,
    intercept = FALSE,
    standardize = FALSE,
    keep = TRUE
  )
  lambda <- cvfit[[args$lambda_rule]]
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("Cross-validation did not return a finite lambda for ", target, ".",
         call. = FALSE)
  }
  coefficients <- .rc_mt_glmnet_coefficients(
    cvfit, lambda, colnames(design)
  )
  effects <- .rc_mt_extract_effective(
    coefficients,
    target_design$edge_metadata$edge_id,
    target_design$condition_levels
  )
  lambda_index <- which.min(abs(cvfit$lambda - lambda))
  prediction <- cvfit$fit.preval[, lambda_index]
  cv_rsq <- .rc_mt_weighted_rsq(y, prediction, weights)

  n_edges <- nrow(target_design$edge_metadata)
  n_conditions <- length(target_design$condition_levels)
  selection_count <- matrix(
    0, nrow = n_edges, ncol = n_conditions,
    dimnames = list(
      target_design$edge_metadata$edge_id,
      target_design$condition_levels
    )
  )
  sign_sum <- selection_count
  successful <- 0L
  if (args$n_stability > 0L) {
    set.seed(target_seed + 1L)
    for (iteration in seq_len(args$n_stability)) {
      selected_rows <- .rc_mt_stability_indices(
        target_design$intercept_stratum, args$stability_fraction
      )
      fit <- tryCatch(
        glmnet::glmnet(
          x = design[selected_rows, , drop = FALSE],
          y = as.numeric(y[selected_rows]),
          family = "gaussian",
          weights = .rc_mt_condition_weights(condition[selected_rows]),
          alpha = args$alpha,
          lambda = lambda,
          penalty.factor = target_design$penalty_factor,
          intercept = FALSE,
          standardize = FALSE
        ),
        error = function(error) NULL
      )
      if (is.null(fit)) next
      one_coefficient <- .rc_mt_glmnet_coefficients(
        fit, lambda, colnames(design)
      )
      one_effect <- .rc_mt_extract_effective(
        one_coefficient,
        target_design$edge_metadata$edge_id,
        target_design$condition_levels
      )$effective
      selected <- abs(one_effect) > args$coefficient_zero_tolerance
      selection_count <- selection_count + selected
      sign_sum <- sign_sum + sign(one_effect) * selected
      successful <- successful + 1L
    }
  }
  if (args$n_stability == 0L) {
    selected <- abs(effects$effective) > args$coefficient_zero_tolerance
    selection_frequency <- selected * 1
    sign_stability <- selected * 1
  } else if (successful > 0L) {
    selection_frequency <- selection_count / successful
    sign_stability <- abs(sign_sum) / pmax(selection_count, 1)
  } else {
    selection_frequency <- selection_count
    sign_stability <- selection_count
  }
  list(
    effects = effects,
    selection_frequency = selection_frequency,
    sign_stability = sign_stability,
    lambda = lambda,
    cv_rsq = cv_rsq,
    n_stability_successful = successful,
    n_stability_requested = args$n_stability,
    edge_metadata = target_design$edge_metadata
  )
}

.rc_mt_format_target_fit <- function(fit, cell_type, condition_col,
                                     celltype_col, args) {
  edge <- fit$edge_metadata
  conditions <- colnames(fit$effects$effective)
  all_rows <- lapply(conditions, function(condition) {
    effective <- fit$effects$effective[, condition]
    frequency <- fit$selection_frequency[, condition]
    stability <- fit$sign_stability[, condition]
    active <- is.finite(fit$cv_rsq) & fit$cv_rsq >= args$min_cv_rsq &
      frequency >= args$min_selection_frequency &
      stability >= args$min_sign_stability &
      abs(effective) >= args$min_abs_effect
    group_values <- data.frame(
      condition_value = condition,
      celltype_value = cell_type,
      stringsAsFactors = FALSE
    )
    names(group_values) <- c(condition_col, celltype_col)
    group_id <- rc_make_stratum_id(group_values, c(condition_col, celltype_col))
    data.frame(
      group_id = rep(group_id, nrow(edge)),
      group_values[rep(1L, nrow(edge)), , drop = FALSE],
      sample_id = rep(group_id, nrow(edge)),
      edge,
      global_estimate = as.numeric(fit$effects$global),
      condition_deviation = as.numeric(fit$effects$deviation[, condition]),
      effective_estimate = as.numeric(effective),
      selection_frequency = as.numeric(frequency),
      sign_stability = as.numeric(stability),
      stability_weight = as.numeric(frequency * stability),
      estimate = as.numeric(effective * frequency * stability),
      active_edge = as.logical(active),
      rsq = rep(fit$cv_rsq, nrow(edge)),
      cv_rsq = rep(fit$cv_rsq, nrow(edge)),
      lambda = rep(fit$lambda, nrow(edge)),
      padj = rep(NA_real_, nrow(edge)),
      evidence_type = rep("multitask_stability_selected", nrow(edge)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  all <- do.call(rbind, all_rows)
  rownames(all) <- NULL
  global <- data.frame(
    cell_type = cell_type,
    edge,
    global_estimate = as.numeric(fit$effects$global),
    cv_rsq = rep(fit$cv_rsq, nrow(edge)),
    lambda = rep(fit$lambda, nrow(edge)),
    stringsAsFactors = FALSE
  )
  names(global)[names(global) == "cell_type"] <- celltype_col
  list(all = all, significant = all[all$active_edge, , drop = FALSE],
       global = global)
}

.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    sample_col = NULL,
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(
      screen_method = "structural",
      min_tf_detection = 0.01,
      min_peak_detection = 0.01,
      min_target_detection = 0.01
    ),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  args <- .rc_mt_validate_args(multitask_args)
  if (!requireNamespace("Pando", quietly = TRUE) ||
      !exists("prepare_grn_design", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "The shared-backbone mode requires Pando >= 1.1.2 from ",
      "1667857557/Pando_regcompass with `prepare_grn_design()`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for multi-task GRN fitting.",
         call. = FALSE)
  }
  pando_install <- .rc_validate_pando_repository()
  if (!is.list(pando_initiate_args) || !is.list(pando_motif_args) ||
      !is.list(pando_design_args)) {
    stop("Pando initiate, motif and design arguments must be lists.",
         call. = FALSE)
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 3 ||
      abs(min_cells - round(min_cells)) > sqrt(.Machine$double.eps)) {
    stop("`min_cells` must be one integer of at least 3.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)
  metadata_columns <- c(condition_col, celltype_col)
  if (!is.null(sample_col)) metadata_columns <- c(metadata_columns, sample_col)
  missing <- setdiff(metadata_columns, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop("Shared-backbone Pando requires cell-type-shared ATAC TF-IDF.",
         call. = FALSE)
  }
  motif_policy <- "user_supplied"
  if (is.null(pfm)) {
    pfm <- .rc_default_pando_motifs()
    motif_policy <- "Pando::motifs"
  }
  if (!length(pfm)) stop("`pfm` must be non-empty.", call. = FALSE)
  region_policy <- "user_supplied"
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
    region_policy <- if (identical(species, "human")) {
      "union(Pando conserved elements, Pando SCREEN ccRE)"
    } else {
      "Pando conserved elements"
    }
  }

  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(
    toupper(.rc_mm_trim_unique(rna_genes)),
    toupper(.rc_mm_trim_unique(metabolic_genes))
  )
  target_genes <- .rc_mm_trim_unique(
    rna_genes[toupper(rna_genes) %in% target_upper]
  )
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"), recursive = TRUE,
             showWarnings = FALSE)
  meta <- object@meta.data
  cell_types <- sort(unique(as.character(meta[[celltype_col]])))
  run_one_celltype <- function(cell_type) {
    cells <- rownames(meta)[as.character(meta[[celltype_col]]) == cell_type]
    condition <- as.character(meta[cells, condition_col, drop = TRUE])
    condition_counts <- table(condition)
    if (any(condition_counts < min_cells)) {
      stop(
        "Cell type ", cell_type, " has conditions below `min_cells`: ",
        paste(names(condition_counts)[condition_counts < min_cells],
              condition_counts[condition_counts < min_cells], sep = "=",
              collapse = "; "),
        call. = FALSE
      )
    }
    obj <- subset(object, cells = cells)
    filtered <- .rc_drop_zero_count_atac_features(
      obj, atac_assay, paste0("Pando cell type ", cell_type)
    )
    obj <- filtered$object
    init_defaults <- list(
      object = obj, peak_assay = atac_assay, rna_assay = rna_assay
    )
    init_defaults[names(pando_initiate_args)] <- NULL
    grn <- do.call(Pando::initiate_grn, c(init_defaults, pando_initiate_args))
    motif_defaults <- list(object = grn, pfm = pfm, genome = genome)
    motif_defaults[names(pando_motif_args)] <- NULL
    grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))
    design_defaults <- list(
      object = grn,
      genes = target_genes,
      screen_group_col = condition_col
    )
    design_defaults[names(pando_design_args)] <- NULL
    design <- do.call(
      Pando::prepare_grn_design, c(design_defaults, pando_design_args)
    )
    Pando::validate_grn_design(design)

    rna <- .rc_pando_assay_data(obj, rna_assay)
    atac <- .rc_pando_assay_data(obj, atac_assay)
    design_cells <- as.character(design$feature_contract$cell_ids)
    if (!identical(design_cells, colnames(rna)) ||
        !identical(design_cells, colnames(atac))) {
      stop("Pando design and normalized assay cells are not aligned.",
           call. = FALSE)
    }
    condition <- as.character(obj@meta.data[design_cells, condition_col,
                                              drop = TRUE])
    sample <- if (is.null(sample_col)) NULL else as.character(
      obj@meta.data[design_cells, sample_col, drop = TRUE]
    )
    candidate <- design$candidate_edges
    target_rows <- split(seq_len(nrow(candidate)), candidate$target)
    target_results <- lapply(names(target_rows), function(target) {
      edge <- candidate[target_rows[[target]], , drop = FALSE]
      target_design <- .rc_mt_target_design(
        rna = rna,
        atac = atac,
        cells = design_cells,
        condition = condition,
        sample = sample,
        edges = edge,
        args = args
      )
      fit <- tryCatch(
        .rc_mt_fit_target(
          y = as.numeric(rna[target, design_cells]),
          target_design = target_design,
          condition = condition,
          args = args,
          target = target
        ),
        error = function(error) error
      )
      if (inherits(fit, "error")) {
        return(list(
          status = data.frame(
            target = target,
            status = "failed",
            error_message = conditionMessage(fit),
            stringsAsFactors = FALSE
          )
        ))
      }
      formatted <- .rc_mt_format_target_fit(
        fit, cell_type, condition_col, celltype_col, args
      )
      list(
        status = data.frame(
          target = target,
          status = "ok",
          error_message = NA_character_,
          n_candidate_edges = nrow(edge),
          n_estimable_edges = nrow(fit$edge_metadata),
          cv_rsq = fit$cv_rsq,
          lambda = fit$lambda,
          n_stability_requested = fit$n_stability_requested,
          n_stability_successful = fit$n_stability_successful,
          stringsAsFactors = FALSE
        ),
        all = formatted$all,
        significant = formatted$significant,
        global = formatted$global,
        edge_metadata = fit$edge_metadata
      )
    })
    status <- do.call(rbind, lapply(target_results, `[[`, "status"))
    failed <- status$status != "ok"
    if (any(failed)) {
      stop(
        "Multi-task target fits failed for cell type ", cell_type, ": ",
        paste(status$target[failed], status$error_message[failed], sep = "=",
              collapse = "; "),
        call. = FALSE
      )
    }
    bind <- function(name) {
      value <- do.call(rbind, lapply(target_results, `[[`, name))
      if (is.null(value)) data.frame() else value
    }
    all <- bind("all")
    significant <- bind("significant")
    global <- bind("global")
    edge_metadata <- bind("edge_metadata")
    candidate <- merge(
      candidate,
      unique(edge_metadata[, c(
        "edge_id", "interaction_scale", "raw_interaction_scale",
        "tf_reference", "estimable"
      ), drop = FALSE]),
      by = "edge_id", all.x = TRUE, sort = FALSE
    )
    candidate[[celltype_col]] <- cell_type
    candidate$design_id <- design$design_id
    status[[celltype_col]] <- cell_type
    if (isTRUE(save_pando_objects)) {
      key <- gsub("[^A-Za-z0-9_.-]+", "_", cell_type)
      saveRDS(grn, file.path(outdir, "pando_objects", paste0(key, ".rds")))
      saveRDS(design, file.path(
        outdir, "pando_objects", paste0(key, "_design.rds")
      ))
    }
    list(
      target_status = status,
      all = all,
      significant = significant,
      global = global,
      candidates = candidate,
      design = design,
      peak_diagnostics = filtered$diagnostics,
      condition_counts = condition_counts
    )
  }

  results <- rc_parallel_lapply(cell_types, run_one_celltype, BPPARAM = BPPARAM)
  bind_results <- function(name) {
    value <- do.call(rbind, lapply(results, `[[`, name))
    if (is.null(value)) data.frame() else value
  }
  all_edges <- bind_results("all")
  significant <- bind_results("significant")
  global <- bind_results("global")
  candidates <- bind_results("candidates")
  target_status <- bind_results("target_status")
  if (!nrow(all_edges)) stop("No multi-task GRN edges were fitted.", call. = FALSE)

  status_rows <- list()
  status_index <- 0L
  for (i in seq_along(results)) {
    cell_type <- cell_types[[i]]
    counts <- results[[i]]$condition_counts
    for (condition in names(counts)) {
      status_index <- status_index + 1L
      values <- data.frame(
        condition_value = condition,
        celltype_value = cell_type,
        stringsAsFactors = FALSE
      )
      names(values) <- c(condition_col, celltype_col)
      group_id <- rc_make_stratum_id(values, c(condition_col, celltype_col))
      status_rows[[status_index]] <- data.frame(
        group_id = group_id,
        values,
        n_cells = as.integer(counts[[condition]]),
        n_target_genes = length(target_genes),
        n_edges = sum(all_edges$group_id == group_id),
        n_significant_edges = sum(significant$group_id == group_id),
        status = "ok",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  sample_status <- do.call(rbind, status_rows)
  rownames(sample_status) <- NULL

  .rc_mm_write_tsv_gz(
    sample_status, file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    target_status, file.path(outdir, "multitask_target_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    candidates, file.path(outdir, "pando_tf_peak_gene_candidates.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    global, file.path(outdir, "pando_tf_peak_gene_global.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    all_edges, file.path(outdir, "pando_tf_peak_gene_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    significant, file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz")
  )

  answer <- list(
    schema_version = "regcompass_multitask_grn_v1",
    grn_mode = "multitask_shared_backbone",
    pando_installed_version = pando_install$version,
    pando_installation = pando_install,
    target_metabolic_genes = target_genes,
    sample_status = sample_status,
    celltype_fit_status = target_status,
    tf_peak_gene_candidates = candidates,
    tf_peak_gene_global = global,
    tf_peak_gene_condition_all = all_edges,
    tf_peak_gene_all = all_edges,
    tf_peak_gene_significant = significant,
    group_cols = c(condition_col, celltype_col),
    design_ids = stats::setNames(
      vapply(results, function(x) x$design$design_id, character(1)),
      cell_types
    ),
    normalization_policy = list(
      rna = "global single-cell NormalizeData",
      atac = "cell-type-shared TF-IDF across conditions",
      candidate_universe = "one Pando structural design per cell type",
      condition_balance = "equal total loss weight per condition",
      pando_motifs = motif_policy,
      pando_regions = region_policy
    ),
    multitask_model = list(
      effective_formula = "theta_edge_condition = global_edge + condition_edge",
      canonical_global = "unweighted mean of effective condition coefficients",
      canonical_deviation = "effective condition coefficient minus canonical global",
      condition_intercepts = if (is.null(sample_col)) {
        "unpenalized condition intercepts"
      } else {
        "unpenalized condition-by-sample intercepts"
      },
      interaction = "normalized TF RNA multiplied by normalized peak ATAC",
      interaction_scale = "root mean within-condition variance",
      stability = "stratified subsampling at fixed cross-validated lambda",
      args = args
    )
  )
  class(answer) <- c("regcompass_multitask_grn", "list")
  saveRDS(answer, file.path(outdir, "multitask_grn_result.rds"))
  answer
}
