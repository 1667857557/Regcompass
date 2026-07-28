# Group-based legacy Pando extraction loaded after pando_grn.R.

rc_extract_pando_tf_peak_gene <- function(
    grn_object,
    group_id = NULL,
    ...,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE) {
  aliases <- list(...)
  if (is.null(group_id) && length(aliases) == 1L) {
    group_id <- aliases[[1L]]
  } else if (length(aliases)) {
    stop("Only one group identifier may be supplied.", call. = FALSE)
  }
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj
  )
  if (!is.character(group_id) || length(group_id) != 1L ||
      is.na(group_id) || !nzchar(trimws(group_id))) {
    stop("`group_id` must be one non-empty character value.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  coefs <- as.data.frame(stats::coef(grn_object), stringsAsFactors = FALSE)
  if (!nrow(coefs)) {
    empty <- data.frame(
      group_id = character(),
      tf = character(),
      target = character(),
      region = character(),
      stringsAsFactors = FALSE
    )
    return(list(all = empty, significant = empty))
  }
  required <- c("tf", "target", "region", "estimate")
  missing <- setdiff(required, colnames(coefs))
  if (length(missing)) {
    stop(
      "Pando coefficient table is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  fit <- tryCatch(
    as.data.frame(Pando::gof(grn_object), stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  if (!nrow(fit) || !all(c("target", "rsq") %in% colnames(fit))) {
    stop("Pando target-model GOF must contain `target` and `rsq`.",
         call. = FALSE)
  }
  keep_fit <- setdiff(
    colnames(fit),
    intersect(colnames(fit), setdiff(colnames(coefs), "target"))
  )
  coefs <- merge(
    coefs,
    fit[, keep_fit, drop = FALSE],
    by = "target",
    all.x = TRUE,
    sort = FALSE
  )
  coefs$group_id <- as.character(group_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[
    , c("group_id", setdiff(colnames(coefs), "group_id")), drop = FALSE
  ]

  estimate <- suppressWarnings(as.numeric(coefs$estimate))
  rsq <- suppressWarnings(as.numeric(coefs$rsq))
  keep <- is.finite(estimate) & abs(estimate) >= min_abs_estimate &
    is.finite(rsq) & rsq >= min_model_rsq
  if ("padj" %in% colnames(coefs)) {
    padj <- suppressWarnings(as.numeric(coefs$padj))
    keep <- keep & is.finite(padj) & padj <= padj_threshold
  } else if (isTRUE(require_padj)) {
    stop(
      paste(
        "Pando network does not contain `padj`; use a p-value-producing",
        "model such as `method = 'glm'`, or set `require_padj = FALSE`."
      ),
      call. = FALSE
    )
  }
  list(all = coefs, significant = coefs[keep, , drop = FALSE])
}

.rc_run_condition_single_cell_grns_legacy <- .rc_run_condition_single_cell_grns

.rc_run_condition_single_cell_grns <- function(...) {
  answer <- .rc_run_condition_single_cell_grns_legacy(...)
  if (is.data.frame(answer$sample_status)) {
    answer$group_status <- answer$sample_status
    answer$sample_status <- NULL
  }
  answer
}
