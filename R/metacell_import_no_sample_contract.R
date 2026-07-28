# Final metacell-import contract. The historical sample-aware builder remains
# outside the active public path, but R CMD check still parses its call site.
# Accepting `...` here prevents a false unused-argument diagnostic while this
# wrapper explicitly rejects sample/pool semantics at runtime.

.rc_import_supercell2_metacells_current_core <-
  rc_import_supercell2_metacells

rc_import_supercell2_metacells <- function(
    metacell_dirs, rna_assay = "RNA", atac_assay = "ATAC",
    condition_col = "condition", celltype_col = "cell_type",
    require_fragments = FALSE, ...) {
  dots <- list(...)
  named <- names(dots) %||% rep("", length(dots))
  forbidden <- named[grepl("sample|pool", named, ignore.case = TRUE)]
  if (length(forbidden)) {
    stop(
      paste(
        "Historical sample/pool arguments are not supported by the current",
        "RegCompass metacell import contract:"
      ),
      paste(unique(forbidden), collapse = ", "),
      call. = FALSE
    )
  }
  if (length(dots)) {
    stop(
      "Unsupported metacell import arguments: ",
      paste(named[nzchar(named)], collapse = ", "),
      call. = FALSE
    )
  }
  .rc_import_supercell2_metacells_current_core(
    metacell_dirs = metacell_dirs,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    condition_col = condition_col,
    celltype_col = celltype_col,
    require_fragments = require_fragments
  )
}
