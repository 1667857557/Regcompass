# RegCompass requires a categorical label because downstream GRN, GPR, and
# reaction comparisons are cell-type resolved. The underlying SuperCell2 label
# argument remains the only grouping information passed into SuperCell2 itself.

.rc_make_supercell2_metacells_unlabelled_core <- rc_make_supercell2_metacells

#' Build label-preserving SuperCell2 metacells for RegCompass
#'
#' @inheritParams rc_make_supercell2_metacells
#' @export
rc_make_supercell2_metacells <- function(
    object, outdir, strata_cols = "condition", label_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    rna_reduction = "pca", atac_reduction = "lsi",
    rna_dims = 1:30, atac_dims = 2:30,
    gamma = 30, seed = 12345L,
    min_cells_per_stratum = 100L, min_metacell_size = 20L,
    min_metacells_per_stratum = 2L,
    fragment_files = FALSE, bgzip_path = "bgzip", tabix_path = "tabix",
    fragment_nb_cl = 1L, save_metacell_object = TRUE, save_counts = TRUE,
    overwrite = FALSE, BPPARAM = NULL,
    on_stratum_error = c("record", "stop"),
    call_peaks_from_fragments = TRUE, macs2_path = NULL,
    peak_calling_effective_genome_size = NULL,
    peak_calling_args = list(), ...) {
  if (!is.character(label_col) || length(label_col) != 1L ||
      is.na(label_col) || !nzchar(trimws(label_col))) {
    stop(
      paste(
        "`label_col` must name one complete categorical metadata column.",
        "RegCompass performs cell-type-resolved GRN and metabolic analysis and",
        "therefore does not support unlabeled metacells."
      ),
      call. = FALSE
    )
  }
  if (!is.logical(call_peaks_from_fragments) ||
      length(call_peaks_from_fragments) != 1L ||
      is.na(call_peaks_from_fragments)) {
    stop("`call_peaks_from_fragments` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(peak_calling_args)) {
    stop("`peak_calling_args` must be a list.", call. = FALSE)
  }
  fragment_enabled <- !identical(fragment_files, FALSE) && !is.null(fragment_files)
  if (fragment_enabled && !isTRUE(save_metacell_object)) {
    stop(
      "Fragment recounting requires `save_metacell_object = TRUE`.",
      call. = FALSE
    )
  }

  built <- .rc_make_supercell2_metacells_unlabelled_core(
    object = object,
    outdir = outdir,
    strata_cols = strata_cols,
    label_col = trimws(label_col),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    rna_reduction = rna_reduction,
    atac_reduction = atac_reduction,
    rna_dims = rna_dims,
    atac_dims = atac_dims,
    gamma = gamma,
    seed = seed,
    min_cells_per_stratum = min_cells_per_stratum,
    min_metacell_size = min_metacell_size,
    min_metacells_per_stratum = min_metacells_per_stratum,
    fragment_files = fragment_files,
    bgzip_path = bgzip_path,
    tabix_path = tabix_path,
    fragment_nb_cl = fragment_nb_cl,
    save_metacell_object = save_metacell_object,
    save_counts = save_counts,
    overwrite = overwrite,
    BPPARAM = BPPARAM,
    on_stratum_error = on_stratum_error,
    ...
  )

  if (!fragment_enabled) {
    built$atac_count_source <- "aggregated_object_peak_counts"
    built$atac_peak_source <- "existing_object_peak_ranges"
    return(built)
  }
  if (!is.data.frame(built$fragment_manifest) ||
      !nrow(built$fragment_manifest)) {
    stop(
      "SuperCell2 fragment aggregation completed without a usable manifest.",
      call. = FALSE
    )
  }
  object_files <- as.character(built$metacell_objects)
  if (!length(object_files)) {
    stop("Fragment recounting requires saved metacell objects.", call. = FALSE)
  }

  for (object_file in object_files) {
    mc <- readRDS(object_file)
    stratum_dir <- dirname(object_file)
    manifest_i <- built$fragment_manifest
    if ("stratum_dir" %in% colnames(manifest_i)) {
      manifest_i <- manifest_i[
        normalizePath(manifest_i$stratum_dir, mustWork = FALSE) ==
          normalizePath(stratum_dir, mustWork = FALSE),
        , drop = FALSE
      ]
    } else {
      manifest_i <- manifest_i[
        as.character(manifest_i$object_cell) %in% colnames(mc),
        , drop = FALSE
      ]
    }
    manifest_i <- .rc_expand_fragment_manifest(manifest_i, colnames(mc))
    expected_peak_source <- if (isTRUE(call_peaks_from_fragments)) {
      "de_novo_macs2_from_metacell_fragments"
    } else {
      "existing_object_peak_ranges"
    }
    already_recounted <- identical(
      tryCatch(mc@misc$atac_count_source, error = function(e) NULL),
      "recomputed_from_metacell_fragments"
    ) && identical(
      tryCatch(mc@misc$atac_peak_source, error = function(e) NULL),
      expected_peak_source
    )
    if (!already_recounted || isTRUE(overwrite)) {
      mc <- .rc_recount_atac_from_fragment_manifest(
        object = mc,
        fragment_manifest = manifest_i,
        atac_assay = atac_assay,
        require_complete = TRUE,
        call_peaks = call_peaks_from_fragments,
        macs2_path = macs2_path,
        effective_genome_size = peak_calling_effective_genome_size,
        peak_calling_args = peak_calling_args,
        peak_calling_outdir = file.path(stratum_dir, "peaks", "macs2")
      )
      saveRDS(mc, object_file)
      saveRDS(
        .rc_as_sparse(.rc_get_assay_counts_safe(mc, atac_assay)),
        file.path(stratum_dir, "atac_counts.rds")
      )
      .rc_write_tsv_gz(
        manifest_i,
        file.path(stratum_dir, "fragments", "fragment_manifest.tsv.gz")
      )
      if (identical(
        tryCatch(mc@misc$atac_peak_source, error = function(e) NULL),
        "de_novo_macs2_from_metacell_fragments"
      )) {
        peak_dir <- file.path(stratum_dir, "peaks")
        dir.create(peak_dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(
          methods::slot(mc[[atac_assay]], "ranges"),
          file.path(peak_dir, "called_peaks.rds")
        )
      }
    }
  }

  refreshed <- rc_import_supercell2_metacells(
    unique(dirname(object_files)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    condition_col = as.character(strata_cols)[[1L]],
    celltype_col = trimws(label_col),
    require_fragments = TRUE
  )
  refreshed$stratum_status <- built$stratum_status
  refreshed$strata_cols <- as.character(strata_cols)
  refreshed$label_col <- trimws(label_col)
  refreshed$construction_contract <-
    "regcompass_split_strata_then_supercell2_label_fragment_recount"
  refreshed$atac_count_source <- "recomputed_from_metacell_fragments"
  refreshed$atac_peak_source <- if (isTRUE(call_peaks_from_fragments)) {
    "de_novo_macs2_from_metacell_fragments"
  } else {
    "existing_object_peak_ranges"
  }
  refreshed
}
