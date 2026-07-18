#' Check biological replicate support for differential testing
rc_check_replicate_design <- function(unit_meta, condition_col = "condition", sample_col = "sample_id", min_samples_per_condition = 2L, strict = TRUE) {
  if (!is.data.frame(unit_meta)) stop("`unit_meta` must be a data.frame.", call. = FALSE)
  missing <- setdiff(c(condition_col, sample_col), colnames(unit_meta))
  if (length(missing) > 0L) stop("`unit_meta` is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  x <- unique(unit_meta[, c(condition_col, sample_col), drop = FALSE])
  x <- x[!is.na(x[[condition_col]]) & !is.na(x[[sample_col]]), , drop = FALSE]
  tab <- table(x[[condition_col]])
  ok <- length(tab) >= 2L && all(tab >= as.integer(min_samples_per_condition))
  if (!ok) {
    msg <- paste0("Insufficient biological replication for differential testing: ", paste(names(tab), as.integer(tab), sep = "=", collapse = ", "), ". Use descriptive summaries only, or set `strict = FALSE` for workflow testing.")
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(tab)
}

.rc_describe_microcompass_by_group_engine <- function(result,
                                             sample_col = "sample_id",
                                             condition_col = "condition",
                                             celltype_col = "cell_type") {
  S <- as.matrix(result$score)
  meta <- result$unit_meta
  meta <- meta[match(colnames(S), as.character(meta$unit_id)), , drop = FALSE]
  parsed <- rc_parse_microcompass_row_id(rownames(S))
  row_meta <- data.frame(row_id = rownames(S), parsed, stringsAsFactors = FALSE)
  groups <- unique(meta[, intersect(c(condition_col, celltype_col), colnames(meta)), drop = FALSE])
  pieces <- lapply(seq_len(nrow(row_meta)), function(i) {
    do.call(rbind, lapply(seq_len(nrow(groups)), function(g) {
      keep <- rep(TRUE, nrow(meta))
      for (nm in colnames(groups)) keep <- keep & as.character(meta[[nm]]) == as.character(groups[[nm]][[g]])
      vals <- as.numeric(S[i, keep])
      data.frame(
        reaction_id = row_meta$reaction_id[i],
        target_direction = row_meta$target_direction[i],
        medium_scenario = row_meta$medium_scenario[i],
        cell_type = if (celltype_col %in% colnames(groups)) as.character(groups[[celltype_col]][[g]]) else NA_character_,
        condition = if (condition_col %in% colnames(groups)) as.character(groups[[condition_col]][[g]]) else NA_character_,
        median_score = stats::median(vals, na.rm = TRUE),
        mean_score = mean(vals, na.rm = TRUE),
        IQR = stats::IQR(vals, na.rm = TRUE),
        n_metacells = sum(keep),
        n_cells = if ("n_cells" %in% colnames(meta)) sum(meta$n_cells[keep], na.rm = TRUE) else NA_real_,
        n_biological_samples = if (sample_col %in% colnames(meta)) length(unique(meta[[sample_col]][keep])) else NA_integer_,
        p_value = NA_real_,
        FDR = NA_real_,
        model_status = "descriptive_only",
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, pieces)
}

#' Test sample-level microCOMPASS differential scores
rc_describe_microcompass_by_group <- function(
    result,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type") {
  metric <- if (!is.null(result$penalty)) "penalty" else "score"
  input <- result
  if (identical(metric, "penalty")) input$score <- result$penalty
  output <- .rc_describe_microcompass_by_group_engine(
    input,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col
  )
  output$analysis_metric <- metric
  output$metric_direction <- if (identical(metric, "penalty")) {
    "lower_is_more_supported"
  } else {
    "higher_is_more_supported"
  }
  output
}

.rc_add_formula_covariates <- function(formula, covariates) {
  formula <- stats::as.formula(formula)
  covariates <- unique(as.character(covariates %||% character()))
  covariates <- covariates[!is.na(covariates) & nzchar(covariates)]
  if (!length(covariates)) return(formula)
  existing <- attr(stats::terms(formula), "term.labels")
  additional <- setdiff(covariates, existing)
  if (!length(additional)) return(formula)
  stats::update.formula(formula, paste(". ~ . +", paste(additional, collapse = " + ")))
}

.rc_aggregate_microcompass_samples <- function(score, meta, sample_col,
                                                condition_col, covariates = NULL) {
  score <- as.matrix(score)
  fields <- unique(c(sample_col, condition_col, covariates))
  missing <- setdiff(fields, colnames(meta))
  if (length(missing)) {
    stop("Metadata is missing sample-level fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  samples <- unique(as.character(meta[[sample_col]]))
  if (anyNA(samples) || any(!nzchar(samples))) {
    stop("Biological sample IDs must be non-missing and non-empty.", call. = FALSE)
  }
  sample_meta <- lapply(samples, function(sample) {
    rows <- which(as.character(meta[[sample_col]]) == sample)
    bad <- vapply(fields, function(field) {
      values <- meta[[field]][rows]
      values <- values[!is.na(values)]
      length(unique(values)) != 1L
    }, logical(1))
    if (any(bad)) {
      stop("Sample `", sample, "` has inconsistent metadata: ",
           paste(fields[bad], collapse = ", "), call. = FALSE)
    }
    meta[rows[[1L]], fields, drop = FALSE]
  })
  sample_meta <- do.call(rbind, sample_meta)
  rownames(sample_meta) <- samples
  aggregated <- do.call(cbind, lapply(samples, function(sample) {
    columns <- which(as.character(meta[[sample_col]]) == sample)
    matrixStats::rowMedians(score[, columns, drop = FALSE], na.rm = TRUE)
  }))
  rownames(aggregated) <- rownames(score)
  colnames(aggregated) <- samples
  list(score = aggregated, meta = sample_meta)
}

.rc_test_microcompass_differential_engine <- function(
    result, formula = score ~ condition,
    method = c("lm", "limma_continuous", "wilcox"),
    sample_col = "sample_id", celltype_col = "cell_type",
    condition_col = "condition", covariates = NULL,
    min_samples_per_group = 3, preferred_min_samples_per_group = 5,
    p_adjust_method = "BH", strict_replicate_design = TRUE,
    test_type = c("omnibus", "pairwise")) {
  method <- match.arg(method)
  test_type <- match.arg(test_type)
  formula <- .rc_add_formula_covariates(formula, covariates)
  if (identical(method, "wilcox") && length(covariates %||% character())) {
    stop("Wilcoxon testing cannot adjust covariates; use lm or limma_continuous.", call. = FALSE)
  }
  score <- as.matrix(result$score)
  meta <- result$unit_meta
  if (is.null(meta) || !"unit_id" %in% colnames(meta)) {
    stop("`result$unit_meta` must contain `unit_id`.", call. = FALSE)
  }
  meta <- meta[match(colnames(score), as.character(meta$unit_id)), , drop = FALSE]
  if (anyNA(meta$unit_id)) stop("`result$unit_meta` is missing rows for score columns.", call. = FALSE)
  required <- c(sample_col, condition_col, covariates)
  missing <- setdiff(required, colnames(meta))
  if (length(missing)) stop("Missing differential metadata: ", paste(missing, collapse = ", "), call. = FALSE)

  parsed <- rc_parse_microcompass_row_id(rownames(score))
  row_meta <- data.frame(row_id = rownames(score), parsed, stringsAsFactors = FALSE)
  celltypes <- if (celltype_col %in% colnames(meta)) unique(as.character(meta[[celltype_col]])) else "all"
  pieces <- lapply(celltypes, function(cell_type) {
    columns <- if (identical(cell_type, "all")) seq_len(ncol(score)) else
      which(as.character(meta[[celltype_col]]) == cell_type)
    if (!length(columns)) return(NULL)
    aggregated <- .rc_aggregate_microcompass_samples(
      score[, columns, drop = FALSE], meta[columns, , drop = FALSE],
      sample_col, condition_col, covariates
    )
    sample_score <- aggregated$score
    sample_meta <- aggregated$meta
    n_by_group <- table(sample_meta[[condition_col]])
    low_power <- length(n_by_group) < 2L || any(n_by_group < min_samples_per_group)
    preferred_low <- length(n_by_group) < 2L || any(n_by_group < preferred_min_samples_per_group)
    if (low_power) {
      if (isTRUE(strict_replicate_design)) {
        stop("Insufficient independent biological samples within cell type `",
             cell_type, "`. Metacells are not biological replicates.", call. = FALSE)
      }
      descriptive_only <- length(n_by_group) < 2L
      return(data.frame(
        reaction_id = row_meta$reaction_id,
        target_direction = row_meta$target_direction,
        cell_type = cell_type,
        medium_scenario = row_meta$medium_scenario,
        contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_,
        p_value = NA_real_,
        n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
        n_biological_samples = nrow(sample_meta),
        method = method, low_sample_power_flag = TRUE,
        preferred_sample_power_flag = preferred_low,
        model_status = if (descriptive_only) "descriptive_only" else "low_sample_power",
        stringsAsFactors = FALSE
      ))
    }
    if (identical(method, "limma_continuous")) {
      return(rc_microcompass_limma(
        sample_score, sample_meta, row_meta, formula,
        sample_col, celltype_col, condition_col,
        min_samples_per_group, preferred_min_samples_per_group,
        p_adjust_method, cell_type, test_type, strict_replicate_design
      ))
    }
    rows <- lapply(seq_len(nrow(sample_score)), function(i) {
      data <- sample_meta
      data$score <- as.numeric(sample_score[i, rownames(sample_meta)])
      fit <- rc_microcompass_fit_one(
        data, formula, method, condition_col, FALSE, test_type = test_type
      )
      data.frame(
        reaction_id = row_meta$reaction_id[[i]],
        target_direction = row_meta$target_direction[[i]],
        cell_type = cell_type,
        medium_scenario = row_meta$medium_scenario[[i]],
        contrast = fit$contrast, effect_size = fit$effect_size,
        statistic = fit$statistic, p_value = fit$p_value,
        n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
        method = method, low_sample_power_flag = FALSE,
        preferred_sample_power_flag = preferred_low,
        model_status = fit$model_status, stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (!length(pieces)) return(data.frame())
  output <- do.call(rbind, pieces)
  adjustment_group <- interaction(
    output$cell_type, output$target_direction,
    output$medium_scenario, output$contrast,
    drop = TRUE, lex.order = TRUE
  )
  output$FDR <- stats::ave(output$p_value, adjustment_group,
                    FUN = function(p) stats::p.adjust(p, method = p_adjust_method))
  rownames(output) <- NULL
  output
}

rc_test_microcompass_differential <- function(
    result, formula = score ~ condition,
    method = c("lm", "limma_continuous", "wilcox"),
    sample_col = "sample_id", celltype_col = "cell_type",
    condition_col = "condition", covariates = NULL,
    min_samples_per_group = 3, preferred_min_samples_per_group = 5,
    p_adjust_method = "BH", strict_replicate_design = TRUE,
    test_type = c("omnibus", "pairwise")) {
  metric <- if (!is.null(result$penalty)) "penalty" else "score"
  input <- result
  if (identical(metric, "penalty")) input$score <- result$penalty
  output <- .rc_test_microcompass_differential_engine(
    input,
    formula = formula,
    method = method,
    sample_col = sample_col,
    celltype_col = celltype_col,
    condition_col = condition_col,
    covariates = covariates,
    min_samples_per_group = min_samples_per_group,
    preferred_min_samples_per_group = preferred_min_samples_per_group,
    p_adjust_method = p_adjust_method,
    strict_replicate_design = strict_replicate_design,
    test_type = test_type
  )
  if (nrow(output)) {
    output$analysis_metric <- metric
    output$metric_direction <- if (identical(metric, "penalty")) {
      "positive_effect_means_higher_penalty_and_weaker_support"
    } else {
      "positive_effect_means_higher_relative_support"
    }
  }
  output
}

rc_microcompass_fit_one <- function(df, formula, method, condition_col,
                                    low_power, test_type = "omnibus") {
  if (isTRUE(low_power)) {
    return(list(contrast = NA_character_, effect_size = NA_real_,
                statistic = NA_real_, p_value = NA_real_,
                model_status = "low_sample_power"))
  }
  if (method == "wilcox") {
    groups <- split(df$score, df[[condition_col]])
    if (length(groups) != 2L) {
      return(list(contrast = NA_character_, effect_size = NA_real_,
                  statistic = NA_real_, p_value = NA_real_,
                  model_status = "requires_two_groups"))
    }
    test <- tryCatch(stats::wilcox.test(groups[[2L]], groups[[1L]]), error = function(e) e)
    if (inherits(test, "error")) {
      return(list(contrast = paste(names(groups), collapse = " vs "),
                  effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_,
                  model_status = conditionMessage(test)))
    }
    return(list(
      contrast = paste(names(groups), collapse = " vs "),
      effect_size = mean(groups[[2L]], na.rm = TRUE) - mean(groups[[1L]], na.rm = TRUE),
      statistic = unname(test$statistic), p_value = test$p.value,
      model_status = "ok"
    ))
  }
  fit <- tryCatch(stats::lm(formula, data = df), error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(contrast = NA_character_, effect_size = NA_real_,
                statistic = NA_real_, p_value = NA_real_,
                model_status = conditionMessage(fit)))
  }
  if (length(unique(df[[condition_col]])) > 2L && identical(test_type, "omnibus")) {
    table <- tryCatch(stats::drop1(fit, test = "F"), error = function(e) e)
    if (inherits(table, "error") || !condition_col %in% rownames(table)) {
      return(list(contrast = paste0(condition_col, "_omnibus"),
                  effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_,
                  model_status = "omnibus_failed"))
    }
    return(list(
      contrast = paste0(condition_col, "_omnibus"), effect_size = NA_real_,
      statistic = unname(table[condition_col, "F value"]),
      p_value = unname(table[condition_col, "Pr(>F)"]), model_status = "ok"
    ))
  }
  coefficients <- summary(fit)$coefficients
  pattern <- paste0("^", make.names(condition_col))
  candidates <- rownames(coefficients)[grepl(pattern, rownames(coefficients))]
  if (!length(candidates)) {
    return(list(contrast = NA_character_, effect_size = NA_real_,
                statistic = NA_real_, p_value = NA_real_,
                model_status = "condition_coefficient_missing"))
  }
  index <- candidates[[1L]]
  list(
    contrast = index,
    effect_size = unname(coefficients[index, "Estimate"]),
    statistic = unname(coefficients[index, "t value"]),
    p_value = unname(coefficients[index, "Pr(>|t|)"]),
    model_status = "ok"
  )
}

rc_microcompass_limma <- function(
    Y, meta, row_meta, formula, sample_col, celltype_col, condition_col,
    min_samples_per_group, preferred_min_samples_per_group,
    p_adjust_method, cell_type, test_type = "omnibus",
    strict_replicate_design = TRUE) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(data.frame(
      reaction_id = row_meta$reaction_id,
      target_direction = row_meta$target_direction,
      cell_type = cell_type, medium_scenario = row_meta$medium_scenario,
      contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_,
      p_value = NA_real_, n_samples_per_group = NA_character_,
      method = "limma_continuous", low_sample_power_flag = TRUE,
      preferred_sample_power_flag = TRUE,
      model_status = "limma package not installed", stringsAsFactors = FALSE
    ))
  }
  n_by_group <- table(meta[[condition_col]])
  low_power <- length(n_by_group) < 2L || any(n_by_group < min_samples_per_group)
  if (isTRUE(strict_replicate_design) && low_power) {
    stop("Insufficient independent biological samples within cell type `",
         cell_type, "`.", call. = FALSE)
  }
  design <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data = meta)
  fit <- limma::eBayes(limma::lmFit(Y, design))
  condition_pattern <- paste0("^", make.names(condition_col))
  condition_coefficients <- colnames(design)[grepl(condition_pattern, colnames(design))]
  if (!length(condition_coefficients)) {
    stop("The model matrix contains no coefficient for `condition_col`.", call. = FALSE)
  }
  coefficient <- if (length(condition_coefficients) > 1L && identical(test_type, "omnibus")) {
    condition_coefficients
  } else {
    condition_coefficients[[1L]]
  }
  table <- if (length(coefficient) > 1L) {
    limma::topTableF(fit, coef = coefficient, number = Inf,
                     sort.by = "none", adjust.method = p_adjust_method)
  } else {
    limma::topTable(fit, coef = coefficient, number = Inf,
                    sort.by = "none", adjust.method = p_adjust_method)
  }
  data.frame(
    reaction_id = row_meta$reaction_id,
    target_direction = row_meta$target_direction,
    cell_type = cell_type, medium_scenario = row_meta$medium_scenario,
    contrast = if (length(coefficient) > 1L) paste0(condition_col, "_omnibus") else coefficient,
    effect_size = if ("logFC" %in% colnames(table)) table$logFC else NA_real_,
    statistic = if ("t" %in% colnames(table)) table$t else table$F,
    p_value = table$P.Value,
    n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
    method = "limma_continuous", low_sample_power_flag = low_power,
    preferred_sample_power_flag = any(n_by_group < preferred_min_samples_per_group),
    model_status = if (low_power) "low_sample_power" else "ok",
    stringsAsFactors = FALSE
  )
}
