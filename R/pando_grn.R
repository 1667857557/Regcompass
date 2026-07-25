#' Load the canonical Pando regulatory-region union
#'
#' The canonical RegCompass GRN uses the union of the conserved-element and
#' SCREEN ccRE region sets bundled with the required Pando fork. Users can
#' override this default through `pando_initiate_args$regions`.
.rc_default_pando_regions <- function() {
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  data_environment <- new.env(parent = emptyenv())
  region_names <- c(
    "phastConsElements20Mammals.UCSC.hg38",
    "SCREEN.ccRE.UCSC.hg38"
  )
  utils::data(
    list = region_names,
    package = "Pando",
    envir = data_environment
  )
  missing <- region_names[!vapply(
    region_names,
    exists,
    logical(1),
    envir = data_environment,
    inherits = FALSE
  )]
  if (length(missing)) {
    stop(
      "Installed Pando does not provide required regulatory-region data: ",
      paste(missing, collapse = ", "),
      ". Install 1667857557/Pando_regcompass.",
      call. = FALSE
    )
  }
  phast_cons <- get(region_names[[1L]], envir = data_environment, inherits = FALSE)
  screen_ccre <- get(region_names[[2L]], envir = data_environment, inherits = FALSE)
  if (!inherits(phast_cons, "GenomicRanges") ||
      !inherits(screen_ccre, "GenomicRanges")) {
    stop("Pando regulatory-region data must be GRanges objects.", call. = FALSE)
  }
  BiocGenerics::union(phast_cons, screen_ccre)
}

#' Extract and filter a Pando TF-peak-gene coefficient table
rc_extract_pando_tf_peak_gene <- function(
    grn_object,
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
      sample_id = character(),
      tf = character(),
      target = character(),
      region = character(),
      stringsAsFactors = FALSE
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
    error = function(error) data.frame()
  )
  if (nrow(fit) && "target" %in% colnames(fit)) {
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
  }
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[, c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE]

  keep <- rep(TRUE, nrow(coefs))
  if ("estimate" %in% colnames(coefs)) {
    estimate <- suppressWarnings(as.numeric(coefs$estimate))
    keep <- keep & is.finite(estimate) & abs(estimate) >= min_abs_estimate
  }
  if ("rsq" %in% colnames(coefs)) {
    rsq <- suppressWarnings(as.numeric(coefs$rsq))
    keep <- keep & is.finite(rsq) & rsq >= min_model_rsq
  } else {
    keep <- rep(FALSE, nrow(coefs))
  }
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
