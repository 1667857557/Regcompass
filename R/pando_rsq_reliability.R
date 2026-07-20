.rc_pando_rsq_is_reliable <- function(rsq, min_model_rsq) {
  if (!is.numeric(min_model_rsq) || length(min_model_rsq) != 1L ||
      !is.finite(min_model_rsq)) {
    stop("`min_model_rsq` must be one finite numeric value.", call. = FALSE)
  }
  value <- suppressWarnings(as.numeric(rsq))
  is.finite(value) & value >= min_model_rsq
}

# Canonical v1.7 correction: missing Pando goodness-of-fit is not evidence.
# The raw coefficient table keeps rsq as NA, while only finite rsq values that
# meet the configured threshold can enter GRN meta-modules or ATAC modifiers.
rc_extract_pando_tf_peak_gene <- function(grn_object,
                                          sample_id,
                                          padj_threshold = 0.05,
                                          min_abs_estimate = 0,
                                          min_model_rsq = 0.1,
                                          require_padj = TRUE) {
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  coefs <- as.data.frame(stats::coef(grn_object), stringsAsFactors = FALSE)
  if (!nrow(coefs)) {
    empty <- data.frame(
      sample_id = character(), tf = character(), target = character(),
      region = character(), stringsAsFactors = FALSE
    )
    return(list(all = empty, significant = empty))
  }
  required <- c("tf", "target", "region")
  missing <- setdiff(required, colnames(coefs))
  if (length(missing)) {
    stop(
      "Pando coefficient table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  fit <- tryCatch(
    as.data.frame(Pando::gof(grn_object), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(fit) && "target" %in% colnames(fit)) {
    keep_fit <- setdiff(
      colnames(fit),
      intersect(colnames(fit), setdiff(colnames(coefs), "target"))
    )
    coefs <- merge(
      coefs, fit[, keep_fit, drop = FALSE],
      by = "target", all.x = TRUE, sort = FALSE
    )
  }
  if (!"rsq" %in% colnames(coefs)) coefs$rsq <- NA_real_
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[, c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE]

  keep <- rep(TRUE, nrow(coefs))
  if ("estimate" %in% colnames(coefs)) {
    keep <- keep & is.finite(coefs$estimate) &
      abs(coefs$estimate) >= min_abs_estimate
  }
  keep <- keep & .rc_pando_rsq_is_reliable(coefs$rsq, min_model_rsq)
  if ("padj" %in% colnames(coefs)) {
    keep <- keep & !is.na(coefs$padj) & coefs$padj <= padj_threshold
  } else if (isTRUE(require_padj)) {
    stop(
      paste0(
        "Pando network does not contain `padj`; use a p-value-producing ",
        "model such as `method = 'glm'`, or set `require_padj = FALSE`."
      ),
      call. = FALSE
    )
  }
  list(all = coefs, significant = coefs[keep, , drop = FALSE])
}
