# Load one Pando dataset through the installed package data index. Never force
# namespace LazyData promises: those promises retain the path of the package
# installation that created them and become invalid after an in-session reinstall.
.rc_pando_data_object <- function(name) {
  if (!is.character(name) || length(name) != 1L ||
      is.na(name) || !nzchar(name)) {
    stop("`name` must be one non-empty character value.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  package_path <- tryCatch(
    getNamespaceInfo(asNamespace("Pando"), "path"),
    error = function(error) ""
  )
  if (!is.character(package_path) || length(package_path) != 1L ||
      is.na(package_path) || !nzchar(package_path) ||
      !dir.exists(package_path)) {
    stop("Cannot resolve the active Pando installation path.", call. = FALSE)
  }

  data_environment <- new.env(parent = emptyenv())
  load_error <- tryCatch({
    utils::data(
      list = name,
      package = "Pando",
      lib.loc = dirname(package_path),
      envir = data_environment
    )
    NULL
  }, error = function(error) error)
  if (inherits(load_error, "error")) {
    stop(
      "Cannot read Pando bundled data `", name, "` from ", package_path,
      ": ", conditionMessage(load_error), ". The active Pando installation ",
      "is incomplete or was replaced after this R session started. Restart R ",
      "and reinstall 1667857557/Pando_regcompass >= 2.0.2.",
      call. = FALSE
    )
  }
  if (!exists(name, envir = data_environment, inherits = FALSE)) {
    stop(
      "Installed Pando does not provide required data object `", name,
      "`. Install 1667857557/Pando_regcompass >= 2.0.2.",
      call. = FALSE
    )
  }
  get(name, envir = data_environment, inherits = FALSE)
}

#' Load the canonical Pando motif collection
#'
#' The canonical RegCompass GRN uses the `motifs` data object bundled with the
#' required Pando fork. Users can override this default by supplying `pfm`.
.rc_default_pando_motifs <- function() {
  motifs <- .rc_pando_data_object("motifs")
  if (is.null(motifs) || !length(motifs)) {
    stop("Pando `motifs` must be a non-empty motif collection.", call. = FALSE)
  }
  motifs
}

.rc_default_pando_motif2tf <- function() {
  motif2tf <- .rc_pando_data_object("motif2tf")
  if (!is.data.frame(motif2tf) || ncol(motif2tf) < 2L) {
    stop(
      "Pando `motif2tf` must be a data frame with at least two columns.",
      call. = FALSE
    )
  }
  motif2tf
}

# Keep the RegCompass Pando route memory bounded. Exact motif-hit coordinates
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
  phast_cons <- .rc_pando_data_object(
    "phastConsElements20Mammals.UCSC.hg38"
  )
  if (!methods::is(phast_cons, "GenomicRanges")) {
    stop(
      "Pando phastCons regulatory regions must be a GRanges object.",
      call. = FALSE
    )
  }
  screen_ccre <- .rc_pando_data_object("SCREEN.ccRE.UCSC.hg38")
  if (!methods::is(screen_ccre, "GenomicRanges")) {
    stop(
      "Pando SCREEN ccRE regulatory regions must be a GRanges object.",
      call. = FALSE
    )
  }
  BiocGenerics::union(phast_cons, screen_ccre)
}

# Resolve every package-backed Pando resource before an outer BiocParallel
# dispatch. Workers receive ordinary serialized objects and never dereference
# installation-specific package data promises.
.rc_materialize_pando_resources <- function(
    pfm = NULL, species = c("human", "mouse"),
    pando_initiate_args = list(), pando_motif_args = list()) {
  species <- match.arg(species)
  if (!is.list(pando_initiate_args)) {
    stop("`pando_initiate_args` must be a list.", call. = FALSE)
  }
  if (!is.list(pando_motif_args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  if (is.null(pfm)) pfm <- .rc_default_pando_motifs()
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
  }
  if (!"motif_tfs" %in% names(pando_motif_args) ||
      is.null(pando_motif_args$motif_tfs)) {
    pando_motif_args$motif_tfs <- .rc_default_pando_motif2tf()
  }
  list(
    pfm = pfm,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args
  )
}

.rc_validate_pando_evidence_filters <- function(
    padj_threshold, require_padj) {
  if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
      !is.finite(padj_threshold) || padj_threshold < 0 ||
      padj_threshold > 1) {
    stop("`padj_threshold` must be one finite number in [0, 1].",
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
    grn_object, sample_id, padj_threshold = 0.05,
    require_padj = TRUE) {
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold, require_padj = require_padj
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
      sample_id = character(), tf = character(), target = character(),
      region = character(), stringsAsFactors = FALSE
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
  if (nrow(fit) && "target" %in% colnames(fit)) {
    keep_fit <- setdiff(
      colnames(fit),
      intersect(colnames(fit), setdiff(colnames(coefs), "target"))
    )
    coefs <- merge(
      coefs, fit[, keep_fit, drop = FALSE], by = "target",
      all.x = TRUE, sort = FALSE
    )
  }
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[
    , c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE
  ]
  estimate <- suppressWarnings(as.numeric(coefs$estimate))
  keep <- is.finite(estimate)
  if ("padj" %in% colnames(coefs)) {
    padj <- suppressWarnings(as.numeric(coefs$padj))
    keep <- keep & is.finite(padj) & padj < padj_threshold
  } else if (isTRUE(require_padj)) {
    stop(
      "Pando network does not contain `padj`; use a p-value-producing model.",
      call. = FALSE
    )
  }
  list(all = coefs, significant = coefs[keep, , drop = FALSE])
}
