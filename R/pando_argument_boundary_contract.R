# Protect workflow-controlled Pando inputs after the grounded parameter policy.
# Users may configure biological region and motif options, but cannot replace the
# validated object, assay roles, motif collection/genome, or GEM target genes
# from inside nested argument lists.

.rc_run_celltype_multitask_grns_parameter_policy_core <-
  .rc_run_celltype_multitask_grns

.rc_assert_pando_nested_boundary <- function(bundle, forbidden, label) {
  if (!is.list(bundle)) {
    stop("`", label, "` must be a list.", call. = FALSE)
  }
  supplied <- unique(names(bundle) %||% character())
  supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
  blocked <- intersect(supplied, forbidden)
  if (length(blocked)) {
    stop(
      "`", label, "` cannot override workflow-controlled fields: ",
      paste(blocked, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 100L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_celltype_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  .rc_assert_pando_nested_boundary(
    pando_initiate_args,
    c("object", "peak_assay", "rna_assay"),
    "pando_initiate_args"
  )
  .rc_assert_pando_nested_boundary(
    pando_motif_args,
    c("object", "pfm", "genome", "store_motif_positions"),
    "pando_motif_args"
  )
  .rc_assert_pando_nested_boundary(
    pando_design_args,
    c("object", "genes"),
    "pando_design_args"
  )
  .rc_run_celltype_multitask_grns_parameter_policy_core(
    object = object,
    gem = gem,
    outdir = outdir,
    pfm = pfm,
    genome = genome,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_cells = min_cells,
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_design_args = pando_design_args,
    multitask_args = multitask_args,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    on_celltype_error = on_celltype_error,
    species = species
  )
}
