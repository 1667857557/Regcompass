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
    on_stratum_error = c("record", "stop"), ...) {
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
  .rc_make_supercell2_metacells_unlabelled_core(
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
}
