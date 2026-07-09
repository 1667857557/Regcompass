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
                                             p_adjust_method = "BH") {
  method <- match.arg(method)
  S <- as.matrix(result$score)
  meta <- result$unit_meta
  if (is.null(meta) || !"unit_id" %in% colnames(meta)) stop("`result$unit_meta` must contain `unit_id`.", call. = FALSE)
  meta <- meta[match(colnames(S), as.character(meta$unit_id)), , drop = FALSE]
  if (anyNA(meta$unit_id)) stop("`result$unit_meta` is missing rows for score matrix columns.", call. = FALSE)
  if (!condition_col %in% colnames(meta)) stop("`condition_col` is missing from `result$unit_meta`.", call. = FALSE)
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
      n_by_group <- table(df[[condition_col]])
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
  n_by_group <- table(meta[[condition_col]])
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
