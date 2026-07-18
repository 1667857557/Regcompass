.rc_infer_macs_effective_genome_size <- function(
    object, atac_assay = "ATAC") {
  ranges <- methods::slot(object[[atac_assay]], "ranges")
  genome_labels <- unique(as.character(GenomeInfoDb::genome(ranges)))
  genome_labels <- genome_labels[
    !is.na(genome_labels) & nzchar(genome_labels)
  ]
  label <- paste(genome_labels, collapse = ";")
  if (grepl("(^|[^A-Za-z])(mm[0-9]+|GRCm|mouse)", label,
            ignore.case = TRUE)) {
    return(1.87e9)
  }
  if (grepl("(^|[^A-Za-z])(hg[0-9]+|GRCh|human)", label,
            ignore.case = TRUE)) {
    return(2.7e9)
  }
  warning(
    "Could not infer the MACS effective genome size from the ATAC ranges; ",
    "using 2.7e9. Set `peak_calling_effective_genome_size` explicitly for ",
    "non-human or unannotated genomes.",
    call. = FALSE
  )
  2.7e9
}

.rc_peak_ids <- function(peaks) {
  paste0(
    as.character(GenomeInfoDb::seqnames(peaks)), "-",
    BiocGenerics::start(peaks), "-",
    BiocGenerics::end(peaks)
  )
}

.rc_call_metacell_peaks <- function(
    fragment_files, object, atac_assay = "ATAC",
    call_peaks = TRUE, macs2_path = NULL,
    effective_genome_size = NULL, peak_calling_args = list(),
    peak_calling_outdir = file.path(tempdir(), "regcompass_macs2"),
    call_peaks_fun = NULL) {
  if (!is.logical(call_peaks) || length(call_peaks) != 1L ||
      is.na(call_peaks)) {
    stop("`call_peaks` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(call_peaks)) {
    peaks <- methods::slot(object[[atac_assay]], "ranges")
    counts <- .rc_get_assay_counts_safe(object, atac_assay)
    if (!inherits(peaks, "GRanges") || !length(peaks) ||
        length(peaks) != nrow(counts)) {
      stop(
        "Existing ATAC ranges must be a non-empty GRanges with one range per ",
        "count-matrix row.",
        call. = FALSE
      )
    }
    peak_ids <- .rc_peak_ids(peaks)
    if (anyDuplicated(peak_ids)) {
      stop("Existing ATAC ranges contain duplicated genomic intervals.",
           call. = FALSE)
    }
    names(peaks) <- peak_ids
    return(list(
      peaks = peaks,
      peak_source = "existing_object_peak_ranges"
    ))
  }
  if (!requireNamespace("Signac", quietly = TRUE) &&
      is.null(call_peaks_fun)) {
    stop("Package 'Signac' is required for metacell peak calling.",
         call. = FALSE)
  }
  if (!is.list(peak_calling_args)) {
    stop("`peak_calling_args` must be a list.", call. = FALSE)
  }
  reserved_args <- intersect(
    names(peak_calling_args),
    c("object", "macs2.path", "outdir", "effective.genome.size")
  )
  if (length(reserved_args)) {
    stop(
      "`peak_calling_args` cannot override explicit peak-calling inputs: ",
      paste(reserved_args, collapse = ", "),
      call. = FALSE
    )
  }
  using_signac_call_peaks <- is.null(call_peaks_fun)
  if (using_signac_call_peaks) {
    call_peaks_fun <- getExportedValue("Signac", "CallPeaks")
  }
  if (!is.null(macs2_path) &&
      (length(macs2_path) != 1L || is.na(macs2_path))) {
    stop("`macs2_path` must be NULL or one non-missing path.", call. = FALSE)
  }
  if (is.null(macs2_path) || !nzchar(as.character(macs2_path))) {
    macs2_path <- unname(Sys.which("macs2"))
  }
  if (!nzchar(as.character(macs2_path)) && using_signac_call_peaks) {
    stop(
      "MACS2 was not found. Install MACS2 or pass its executable path in ",
      "`metacell_args$macs2_path`.",
      call. = FALSE
    )
  }
  effective_genome_size <- effective_genome_size %||%
    .rc_infer_macs_effective_genome_size(object, atac_assay = atac_assay)
  if (!is.numeric(effective_genome_size) ||
      length(effective_genome_size) != 1L ||
      !is.finite(effective_genome_size) || effective_genome_size <= 0) {
    stop("`effective_genome_size` must be one positive finite number.",
         call. = FALSE)
  }
  outdir <- peak_calling_outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  default_args <- list(
    object = unique(as.character(fragment_files)),
    macs2.path = macs2_path,
    outdir = outdir,
    broad = FALSE,
    format = "BED",
    effective.genome.size = effective_genome_size,
    extsize = 200,
    shift = -100,
    additional.args = NULL,
    name = paste0(
      "RegCompass_",
      abs(sum(utf8ToInt(paste(colnames(object), collapse = "|"))))
    ),
    cleanup = TRUE,
    verbose = FALSE
  )
  default_args[names(peak_calling_args)] <- NULL
  peaks <- do.call(call_peaks_fun, c(default_args, peak_calling_args))
  if (!inherits(peaks, "GRanges") || !length(peaks)) {
    stop("MACS2 did not return a non-empty GRanges peak set.",
         call. = FALSE)
  }
  peaks <- GenomicRanges::reduce(peaks, ignore.strand = TRUE)
  peaks <- GenomeInfoDb::sortSeqlevels(peaks)
  peaks <- sort(peaks)
  names(peaks) <- .rc_peak_ids(peaks)
  list(
    peaks = peaks,
    peak_source = "de_novo_macs2_from_metacell_fragments",
    effective_genome_size = effective_genome_size,
    macs2_path = as.character(macs2_path)
  )
}

.rc_fragment_objects_and_counts <- function(
    object, fragment_manifest, atac_assay = "ATAC",
    require_complete = TRUE, process_n = 2000L,
    create_fragment_fun = NULL, feature_matrix_fun = NULL,
    call_peaks = TRUE,
    macs2_path = NULL, effective_genome_size = NULL,
    peak_calling_args = list(), peak_calling_outdir = NULL,
    call_peaks_fun = NULL) {
  if (!is.numeric(process_n) || length(process_n) != 1L ||
      !is.finite(process_n) || process_n < 1 ||
      abs(process_n - round(process_n)) > sqrt(.Machine$double.eps)) {
    stop("`process_n` must be one positive integer.", call. = FALSE)
  }
  if (!requireNamespace("Signac", quietly = TRUE) &&
      (is.null(create_fragment_fun) || is.null(feature_matrix_fun) ||
       (isTRUE(call_peaks) && is.null(call_peaks_fun)))) {
    stop(
      "Package 'Signac' is required to call peaks and recount ATAC fragments.",
      call. = FALSE
    )
  }
  if (!atac_assay %in% names(object@assays) ||
      !inherits(object[[atac_assay]], "ChromatinAssay")) {
    stop("ATAC assay `", atac_assay,
         "` must be a Signac ChromatinAssay.", call. = FALSE)
  }
  manifest <- .rc_validate_fragment_recount_manifest(
    fragment_manifest,
    object_cells = colnames(object),
    require_complete = require_complete
  )
  if (!nrow(manifest)) {
    return(list(
      counts = NULL, fragments = list(), manifest = manifest,
      peaks = NULL, peak_source = NA_character_
    ))
  }
  if (is.null(create_fragment_fun)) {
    create_fragment_fun <- getExportedValue("Signac", "CreateFragmentObject")
  }
  if (is.null(feature_matrix_fun)) {
    feature_matrix_fun <- getExportedValue("Signac", "FeatureMatrix")
  }
  files <- unique(manifest$fragment_file)
  peak_result <- .rc_call_metacell_peaks(
    fragment_files = files,
    object = object,
    atac_assay = atac_assay,
    call_peaks = call_peaks,
    macs2_path = macs2_path,
    effective_genome_size = effective_genome_size,
    peak_calling_args = peak_calling_args,
    peak_calling_outdir = peak_calling_outdir %||%
      file.path(tempdir(), "regcompass_macs2"),
    call_peaks_fun = call_peaks_fun
  )
  peaks <- peak_result$peaks
  feature_ids <- .rc_peak_ids(peaks)
  names(peaks) <- feature_ids
  cell_ids <- colnames(object)
  fragment_objects <- vector("list", length(files))
  count_matrices <- vector("list", length(files))
  for (i in seq_along(files)) {
    path <- files[[i]]
    rows <- manifest[manifest$fragment_file == path, , drop = FALSE]
    cell_map <- stats::setNames(
      as.character(rows$fragment_barcode),
      as.character(rows$object_cell)
    )
    fragment_objects[[i]] <- create_fragment_fun(
      path = path,
      cells = cell_map,
      validate.fragments = FALSE
    )
    counts_i <- feature_matrix_fun(
      fragments = list(fragment_objects[[i]]),
      features = peaks,
      cells = names(cell_map),
      process_n = as.integer(process_n),
      verbose = FALSE
    )
    count_matrices[[i]] <- .rc_align_peak_count_matrix(
      counts_i,
      feature_ids = feature_ids,
      cell_ids = cell_ids
    )
  }
  counts <- .rc_as_sparse(Reduce(`+`, count_matrices))
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("Fragment-derived ATAC peak counts contain invalid values.",
         call. = FALSE)
  }
  c(
    list(
      counts = counts,
      fragments = fragment_objects,
      manifest = manifest,
      peaks = peaks,
      peak_source = peak_result$peak_source
    ),
    peak_result[setdiff(names(peak_result), c("peaks", "peak_source"))]
  )
}

.rc_recount_atac_from_fragment_manifest <- function(
    object, fragment_manifest, atac_assay = "ATAC",
    require_complete = TRUE, process_n = 2000L,
    create_fragment_fun = NULL, feature_matrix_fun = NULL,
    call_peaks = TRUE,
    macs2_path = NULL, effective_genome_size = NULL,
    peak_calling_args = list(), peak_calling_outdir = NULL,
    call_peaks_fun = NULL) {
  answer <- .rc_fragment_objects_and_counts(
    object = object,
    fragment_manifest = fragment_manifest,
    atac_assay = atac_assay,
    require_complete = require_complete,
    process_n = process_n,
    create_fragment_fun = create_fragment_fun,
    feature_matrix_fun = feature_matrix_fun,
    call_peaks = call_peaks,
    macs2_path = macs2_path,
    effective_genome_size = effective_genome_size,
    peak_calling_args = peak_calling_args,
    peak_calling_outdir = peak_calling_outdir,
    call_peaks_fun = call_peaks_fun
  )
  if (is.null(answer$counts)) return(object)
  old_assay <- object[[atac_assay]]
  annotation <- tryCatch(
    Signac::Annotation(old_assay),
    error = function(e) NULL
  )
  new_assay <- Signac::CreateChromatinAssay(
    counts = answer$counts,
    ranges = answer$peaks,
    annotation = annotation,
    min.cells = 0,
    min.features = 0
  )
  try(SeuratObject::Key(new_assay) <- SeuratObject::Key(old_assay), silent = TRUE)
  object[[atac_assay]] <- new_assay
  object[[paste0("nCount_", atac_assay)]] <-
    as.numeric(Matrix::colSums(answer$counts))
  object[[paste0("nFeature_", atac_assay)]] <-
    as.numeric(Matrix::colSums(answer$counts > 0))
  mapped <- unlist(lapply(
    split(answer$manifest, answer$manifest$fragment_file),
    function(x) unique(as.character(x$object_cell))
  ), use.names = FALSE)
  if (!anyDuplicated(mapped)) {
    fragment_setter <- get("Fragments<-", envir = asNamespace("Signac"))
    object[[atac_assay]] <- fragment_setter(
      object[[atac_assay]],
      value = answer$fragments
    )
    registration <- "registered"
  } else {
    registration <- "not_registered_overlapping_fragment_files"
  }
  object@misc$atac_count_source <- "recomputed_from_metacell_fragments"
  object@misc$atac_peak_source <- answer$peak_source
  object@misc$atac_fragment_recount <- list(
    n_fragment_files = length(unique(answer$manifest$fragment_file)),
    n_peaks = nrow(answer$counts),
    n_metacells = ncol(answer$counts),
    fragment_registration = registration,
    peak_source = answer$peak_source,
    effective_genome_size = answer$effective_genome_size %||% NA_real_,
    macs2_path = answer$macs2_path %||% NA_character_
  )
  object
}
