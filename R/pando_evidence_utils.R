.rc_pando_assay_data <- function(object, assay) {
  value <- tryCatch(
    .rc_get_assay_matrix(object, assay, "data"),
    error = function(e) NULL
  )
  if (is.null(value) || nrow(value) == 0L || ncol(value) == 0L) {
    stop(
      "Pando evidence projection requires one non-empty normalized `data` ",
      "layer per assay.",
      call. = FALSE
    )
  }
  value
}

.rc_pando_region_key <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^([^:]+):(\\d+)-(\\d+)$", "\\1-\\2-\\3", x)
  x <- sub("^([^:]+):(\\d+):(\\d+)$", "\\1-\\2-\\3", x)
  x
}

.rc_pando_region_coordinates <- function(x) {
  keys <- .rc_pando_region_key(x)
  matched <- regexec("^(.*)-(\\d+)-(\\d+)$", keys)
  parts <- regmatches(keys, matched)
  data.frame(
    seqname = vapply(parts, function(value) {
      if (length(value) == 4L) value[[2L]] else NA_character_
    }, character(1)),
    start = suppressWarnings(as.numeric(vapply(parts, function(value) {
      if (length(value) == 4L) value[[3L]] else NA_character_
    }, character(1)))),
    end = suppressWarnings(as.numeric(vapply(parts, function(value) {
      if (length(value) == 4L) value[[4L]] else NA_character_
    }, character(1)))),
    stringsAsFactors = FALSE
  )
}

# Pando reports the regulatory candidate interval in its coefficient table.
# That interval can be the intersection of a supplied regulatory element and an
# assay peak, so it is not necessarily identical to the ATAC feature name. Map
# exact names first and otherwise recover the source peak by genomic overlap.
.rc_map_pando_regions_to_atac <- function(regions, atac_features) {
  region_keys <- .rc_pando_region_key(regions)
  atac_keys <- .rc_pando_region_key(atac_features)
  if (anyDuplicated(atac_keys)) {
    stop("ATAC features are duplicated after Pando region normalization.",
         call. = FALSE)
  }

  answer <- match(region_keys, atac_keys)
  unresolved <- which(is.na(answer))
  if (!length(unresolved)) return(answer)

  candidate <- .rc_pando_region_coordinates(region_keys)
  peaks <- .rc_pando_region_coordinates(atac_keys)
  valid_peaks <- !is.na(peaks$seqname) & is.finite(peaks$start) &
    is.finite(peaks$end) & nzchar(peaks$seqname) &
    peaks$start <= peaks$end
  valid_candidates <- !is.na(candidate$seqname) &
    is.finite(candidate$start) & is.finite(candidate$end) &
    nzchar(candidate$seqname) & candidate$start <= candidate$end
  candidate_index <- unresolved[valid_candidates[unresolved]]
  peak_index <- which(valid_peaks)
  if (!length(candidate_index) || !length(peak_index)) return(answer)

  candidate_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    candidate[candidate_index, , drop = FALSE],
    seqnames.field = "seqname", start.field = "start", end.field = "end"
  )
  peak_ranges <- GenomicRanges::makeGRangesFromDataFrame(
    peaks[peak_index, , drop = FALSE],
    seqnames.field = "seqname", start.field = "start", end.field = "end"
  )
  range_hits <- GenomicRanges::findOverlaps(
    candidate_ranges, peak_ranges, ignore.strand = TRUE
  )
  if (!length(range_hits)) return(answer)
  hit_table <- as.data.frame(range_hits)
  hits_by_candidate <- split(hit_table$subjectHits, hit_table$queryHits)
  for (query in names(hits_by_candidate)) {
    i <- candidate_index[[as.integer(query)]]
    hits <- peak_index[hits_by_candidate[[query]]]
    overlap <- pmin(peaks$end[hits], candidate$end[[i]]) -
      pmax(peaks$start[hits], candidate$start[[i]]) + 1
    # ATAC peak sets should be disjoint. If they are not, prefer the peak with
    # the greatest shared span, then the closest midpoint, then assay order.
    midpoint_distance <- abs(
      (peaks$start[hits] + peaks$end[hits]) -
        (candidate$start[[i]] + candidate$end[[i]])
    )
    answer[[i]] <- hits[order(-overlap, midpoint_distance, hits)[[1L]]]
  }
  answer
}

.rc_pando_projection_from_group_means <- function(
    rna, atac, edges_by_group, cells_by_group, targets) {
  answer <- matrix(
    0, nrow = length(cells_by_group), ncol = length(targets),
    dimnames = list(names(cells_by_group), tolower(targets))
  )
  rna_names <- toupper(colnames(rna))
  atac_features <- colnames(atac)
  for (group in names(cells_by_group)) {
    edge <- edges_by_group[[group]]
    cells <- intersect(cells_by_group[[group]], rownames(rna))
    cells <- intersect(cells, rownames(atac))
    if (is.null(edge) || !nrow(edge) || !length(cells)) next
    tf_index <- match(toupper(as.character(edge$tf)), rna_names)
    region <- if ("atac_feature_id" %in% colnames(edge)) {
      feature_id <- as.character(edge$atac_feature_id)
      fallback <- is.na(feature_id) | !nzchar(trimws(feature_id))
      feature_id[fallback] <- as.character(edge$region[fallback])
      feature_id
    } else {
      as.character(edge$region)
    }
    peak_index <- .rc_map_pando_regions_to_atac(region, atac_features)
    estimate <- suppressWarnings(as.numeric(edge$estimate))
    usable <- !is.na(tf_index) & !is.na(peak_index) & is.finite(estimate)
    if (!any(usable)) next
    edge <- edge[usable, , drop = FALSE]
    tf_index <- tf_index[usable]
    peak_index <- peak_index[usable]
    estimate <- estimate[usable]
    tf_mean <- colMeans(rna[cells, tf_index, drop = FALSE])
    peak_mean <- colMeans(atac[cells, peak_index, drop = FALSE])
    contribution <- estimate * tf_mean * peak_mean
    edge_target <- tolower(as.character(edge$target))
    totals <- tapply(contribution, edge_target, sum)
    selected <- intersect(names(totals), colnames(answer))
    answer[group, selected] <- totals[selected]
  }
  answer
}
