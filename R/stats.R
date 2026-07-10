#' Check biological replicate support for differential testing
#' @export
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

rc_parse_microcompass_row_id <- function(x) {
  core <- sub("::medium=.*$", "", x)
  parts <- strsplit(core, "::", fixed = TRUE)
  medium <- ifelse(grepl("::medium=", x, fixed = TRUE),
                   sub("^.*::medium=", "", x),
                   vapply(parts, function(z) z[[3]] %||% NA_character_, character(1)))
  data.frame(
    reaction_id = vapply(parts, function(z) z[[1]] %||% NA_character_, character(1)),
    target_direction = vapply(parts, function(z) z[[2]] %||% NA_character_, character(1)),
    medium_scenario = medium,
    stringsAsFactors = FALSE
  )
}

rc_describe_microcompass_by_group <- function(result,
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
#' @export
rc_test_microcompass_differential <- function(result,
                                             formula = score ~ condition,
                                             method = c("lm", "limma_continuous", "wilcox"),
                                             sample_col = "sample_id",
                                             celltype_col = "cell_type",
                                             condition_col = "condition",
                                             covariates = NULL,
                                             min_samples_per_group = 3,
                                             preferred_min_samples_per_group = 5,
                                             p_adjust_method = "BH",
                                             strict_replicate_design = TRUE) {
  method <- match.arg(method)
  S <- as.matrix(result$score)
  meta <- result$unit_meta
  if (is.null(meta) || !"unit_id" %in% colnames(meta)) stop("`result$unit_meta` must contain `unit_id`.", call. = FALSE)
  meta <- meta[match(colnames(S), as.character(meta$unit_id)), , drop = FALSE]
  if (anyNA(meta$unit_id)) stop("`result$unit_meta` is missing rows for score matrix columns.", call. = FALSE)
  if (!condition_col %in% colnames(meta)) stop("`condition_col` is missing from `result$unit_meta`.", call. = FALSE)
  if (!sample_col %in% colnames(meta)) stop("`sample_col` is missing from `result$unit_meta`.", call. = FALSE)
  sample_condition <- unique(meta[, c(sample_col, condition_col), drop = FALSE])
  n_samples_by_group <- table(sample_condition[[condition_col]])
  enough_samples <- length(n_samples_by_group) >= 2L && all(n_samples_by_group >= min_samples_per_group)
  if (!enough_samples) {
    if (isTRUE(strict_replicate_design)) {
      stop("Insufficient independent biological samples. Metacells are not biological replicates.", call. = FALSE)
    }
    return(rc_describe_microcompass_by_group(result = result, sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col))
  }
  parsed <- rc_parse_microcompass_row_id(rownames(S))
  row_meta <- data.frame(row_id = rownames(S), parsed, stringsAsFactors = FALSE)
  celltypes <- if (celltype_col %in% colnames(meta)) unique(as.character(meta[[celltype_col]])) else NA_character_
  pieces <- lapply(celltypes, function(ct) {
    cols <- if (!is.na(ct) && celltype_col %in% colnames(meta)) which(as.character(meta[[celltype_col]]) == ct) else seq_len(ncol(S))
    if (length(cols) == 0L) return(NULL)
    if (method == "limma_continuous") {
      return(rc_microcompass_limma(S[, cols, drop = FALSE], meta[cols, , drop = FALSE], row_meta,
                                   formula, sample_col, celltype_col, condition_col,
                                   min_samples_per_group, preferred_min_samples_per_group,
                                   p_adjust_method, ct))
    }
    rows <- lapply(seq_len(nrow(S)), function(i) {
      df <- data.frame(score = as.numeric(S[i, cols]), meta[cols, , drop = FALSE], stringsAsFactors = FALSE)
      sample_condition_i <- unique(df[, c(sample_col, condition_col), drop = FALSE])
      n_by_group <- table(sample_condition_i[[condition_col]])
      low_power <- length(n_by_group) < 2L || any(n_by_group < min_samples_per_group)
      preferred_low <- length(n_by_group) < 2L || any(n_by_group < preferred_min_samples_per_group)
      fit <- rc_microcompass_fit_one(df, formula, method, condition_col, low_power)
      data.frame(
        reaction_id = row_meta$reaction_id[i],
        target_direction = row_meta$target_direction[i],
        cell_type = if (!is.na(ct)) ct else NA_character_,
        medium_scenario = row_meta$medium_scenario[i],
        contrast = fit$contrast,
        effect_size = fit$effect_size,
        statistic = fit$statistic,
        p_value = fit$p_value,
        n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
        method = method,
        low_sample_power_flag = low_power,
        preferred_sample_power_flag = preferred_low,
        model_status = fit$model_status,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  out <- do.call(rbind, pieces)
  if (is.null(out) || nrow(out) == 0L) return(data.frame())
  out$FDR <- stats::p.adjust(out$p_value, method = p_adjust_method)
  out
}

rc_microcompass_fit_one <- function(df, formula, method, condition_col, low_power) {
  if (isTRUE(low_power)) {
    return(list(contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_,
                p_value = NA_real_, model_status = "low_sample_power"))
  }
  if (method == "wilcox") {
    groups <- split(df$score, df[[condition_col]])
    if (length(groups) != 2L) return(list(contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_, model_status = "requires_two_groups"))
    wt <- tryCatch(stats::wilcox.test(groups[[2]], groups[[1]]), error = function(e) e)
    if (inherits(wt, "error")) return(list(contrast = paste(names(groups), collapse = " vs "), effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_, model_status = conditionMessage(wt)))
    return(list(contrast = paste(names(groups), collapse = " vs "),
                effect_size = mean(groups[[2]], na.rm = TRUE) - mean(groups[[1]], na.rm = TRUE),
                statistic = unname(wt$statistic), p_value = wt$p.value, model_status = "ok"))
  }
  fit <- tryCatch(stats::lm(formula, data = df), error = function(e) e)
  if (inherits(fit, "error")) return(list(contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_, model_status = conditionMessage(fit)))
  cf <- summary(fit)$coefficients
  idx <- setdiff(rownames(cf), "(Intercept)")[1]
  if (is.na(idx)) return(list(contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_, p_value = NA_real_, model_status = "no_non_intercept_term"))
  list(contrast = idx, effect_size = unname(cf[idx, "Estimate"]), statistic = unname(cf[idx, "t value"]),
       p_value = unname(cf[idx, "Pr(>|t|)"]), model_status = "ok")
}

rc_microcompass_limma <- function(Y, meta, row_meta, formula, sample_col, celltype_col, condition_col,
                                  min_samples_per_group, preferred_min_samples_per_group,
                                  p_adjust_method, cell_type) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    return(data.frame(reaction_id = row_meta$reaction_id,
                      target_direction = row_meta$target_direction,
                      cell_type = if (!is.na(cell_type)) cell_type else NA_character_,
                      medium_scenario = row_meta$medium_scenario,
                      contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_,
                      p_value = NA_real_, n_samples_per_group = NA_character_,
                      method = "limma_continuous", low_sample_power_flag = TRUE,
                      preferred_sample_power_flag = TRUE,
                      model_status = "limma package not installed",
                      stringsAsFactors = FALSE))
  }
  sample_condition <- unique(meta[, c(sample_col, condition_col), drop = FALSE])
  n_by_group <- table(sample_condition[[condition_col]])
  low_power <- length(n_by_group) < 2L || any(n_by_group < min_samples_per_group)
  preferred_low <- length(n_by_group) < 2L || any(n_by_group < preferred_min_samples_per_group)
  if (low_power) {
    return(data.frame(reaction_id = row_meta$reaction_id,
                      target_direction = row_meta$target_direction,
                      cell_type = if (!is.na(cell_type)) cell_type else NA_character_,
                      medium_scenario = row_meta$medium_scenario,
                      contrast = NA_character_, effect_size = NA_real_, statistic = NA_real_,
                      p_value = NA_real_, n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
                      method = "limma_continuous", low_sample_power_flag = TRUE,
                      preferred_sample_power_flag = preferred_low,
                      model_status = "low_sample_power",
                      stringsAsFactors = FALSE))
  }
  design <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data = meta)
  fit <- limma::eBayes(limma::lmFit(Y, design))
  coef_name <- setdiff(colnames(design), "(Intercept)")[1]
  tab <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none", adjust.method = p_adjust_method)
  data.frame(reaction_id = row_meta$reaction_id,
             target_direction = row_meta$target_direction,
             cell_type = if (!is.na(cell_type)) cell_type else NA_character_,
             medium_scenario = row_meta$medium_scenario,
             contrast = coef_name,
             effect_size = tab$logFC,
             statistic = tab$t,
             p_value = tab$P.Value,
             n_samples_per_group = paste(names(n_by_group), as.integer(n_by_group), sep = "=", collapse = ";"),
             method = "limma_continuous",
             low_sample_power_flag = FALSE,
             preferred_sample_power_flag = preferred_low,
             model_status = "ok",
             stringsAsFactors = FALSE)
}
