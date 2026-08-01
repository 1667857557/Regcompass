#' Load the canonical Pando motif collection
#'
#' The canonical RegCompass GRN uses the `motifs` data object bundled with the
#' required Pando fork. Users can override this default by supplying `pfm`.
.rc_default_pando_motifs <- function() {
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  data_environment <- new.env(parent = emptyenv())
  utils::data(
    list = "motifs",
    package = "Pando",
    envir = data_environment
  )
  if (!exists("motifs", envir = data_environment, inherits = FALSE)) {
    stop(
      paste(
        "Installed Pando does not provide the required `motifs` data object.",
        "Install 1667857557/Pando_regcompass."
      ),
      call. = FALSE
    )
  }
  motifs <- get("motifs", envir = data_environment, inherits = FALSE)
  if (is.null(motifs) || !length(motifs)) {
    stop("Pando `motifs` must be a non-empty motif collection.", call. = FALSE)
  }
  motifs
}

# Keep the RegCompass Pando route memory bounded.  Exact motif-hit coordinates
# are not consumed anywhere downstream; the candidate construction only needs
# the peak-by-motif incidence matrix.
.rc_regcompass_motif_args <- function(args = list()) {
  if (!is.list(args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  aliases <- c(
    "exact_positions", "exact_motif_positions", "keep_motif_positions",
    "store_motif_positions", "return_motif_positions"
  )
  forbidden <- intersect(names(args), aliases)
  if (length(forbidden)) {
    stop(
      "The RegCompass route always skips exact motif positions; remove: ",
      paste(forbidden, collapse = ", "), call. = FALSE
    )
  }
  available <- names(formals(Pando::find_motifs))
  position_arg <- aliases[aliases %in% available][1L]
  if (is.na(position_arg) && "..." %in% available) {
    position_arg <- "exact_positions"
  }
  if (is.na(position_arg)) {
    stop(
      "Installed Pando cannot skip exact motif positions. Install the current ",
      "1667857557/Pando_regcompass release with binary motif storage.",
      call. = FALSE
    )
  }
  args[[position_arg]] <- FALSE
  args
}

#' Load the canonical species-specific Pando regulatory regions
#'
#' Human analyses use the union of the conserved-element and SCREEN ccRE sets.
#' Mouse analyses require user-supplied mouse-coordinate regions through
#' `pando_initiate_args$regions`.
.rc_default_pando_regions <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (identical(species, "mouse")) {
    stop(
      paste(
        "No mouse-coordinate regulatory-region set is bundled.",
        "Supply `pando_initiate_args$regions` in the same mouse genome build",
        "as the ATAC peaks and `genome`."
      ),
      call. = FALSE
    )
  }
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
  phast_cons <- get(
    "phastConsElements20Mammals.UCSC.hg38",
    envir = data_environment,
    inherits = FALSE
  )
  if (!methods::is(phast_cons, "GenomicRanges")) {
    stop(
      "Pando phastCons regulatory regions must be a GRanges object.",
      call. = FALSE
    )
  }
  screen_ccre <- get(
    "SCREEN.ccRE.UCSC.hg38",
    envir = data_environment,
    inherits = FALSE
  )
  if (!methods::is(screen_ccre, "GenomicRanges")) {
    stop(
      "Pando SCREEN ccRE regulatory regions must be a GRanges object.",
      call. = FALSE
    )
  }
  BiocGenerics::union(phast_cons, screen_ccre)
}

.rc_validate_pando_evidence_filters <- function(
    padj_threshold, min_abs_estimate, min_model_rsq, require_padj) {
  if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
      !is.finite(padj_threshold) || padj_threshold < 0 ||
      padj_threshold > 1) {
    stop("`padj_threshold` must be one finite number in [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(min_abs_estimate) || length(min_abs_estimate) != 1L ||
      !is.finite(min_abs_estimate) || min_abs_estimate < 0) {
    stop("`min_abs_estimate` must be one finite non-negative number.",
         call. = FALSE)
  }
  if (!is.numeric(min_model_rsq) || length(min_model_rsq) != 1L ||
      !is.finite(min_model_rsq) || min_model_rsq < 0) {
    stop("`min_model_rsq` must be one finite non-negative number.",
         call. = FALSE)
  }
  if (!is.logical(require_padj) || length(require_padj) != 1L ||
      is.na(require_padj)) {
    stop("`require_padj` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Extract and filter a Pando TF-peak-gene coefficient table
rc_extract_pando_tf_peak_gene <- function(
    grn_object,
    sample_id,
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE) {
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj
  )
  if (!is.character(sample_id) || length(sample_id) != 1L ||
      is.na(sample_id) || !nzchar(trimws(sample_id))) {
    stop("`sample_id` must be one non-empty character value.", call. = FALSE)
  }
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
  required <- c("tf", "target", "region", "estimate")
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
  if (!nrow(fit) || !all(c("target", "rsq") %in% colnames(fit))) {
    stop(
      "Pando target-model GOF must contain `target` and `rsq`.",
      call. = FALSE
    )
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
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[
    , c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE
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
