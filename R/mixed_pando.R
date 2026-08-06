# Cell-type-specific routing between common-dictionary and standard Pando.

.rc_bind_pando_field <- function(results, field) {
  values <- lapply(results, function(x) x[[field]])
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  .rc_bind_frames_fill(values)
}

.rc_merge_pando_results_core <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
  results <- c(
    if (is.null(condition_result)) list() else list(condition_result),
    standard_results
  )
  if (!length(results)) stop("No Pando result was produced.", call. = FALSE)
  standard_objects <- list()
  for (value in standard_results) {
    for (name in names(value$standard_pando_objects)) {
      if (name %in% names(standard_objects)) {
        stop("Duplicated standard-Pando cell type: ", name, call. = FALSE)
      }
      standard_objects[[name]] <- value$standard_pando_objects[[name]]
    }
  }
  condition_fits <- if (is.null(condition_result)) list() else
    condition_result$condition_grn_fits
  paired_meta <- .rc_bind_pando_field(results, "paired_cell_metadata")
  if (nrow(paired_meta)) {
    paired_meta <- paired_meta[!duplicated(paired_meta$cell_id), , drop = FALSE]
  }
  routing <- rbind(
    if (length(condition_types)) data.frame(
      cell_type = condition_types, analysis_mode = "condition_grn",
      stringsAsFactors = FALSE
    ),
    if (length(standard_types)) data.frame(
      cell_type = standard_types, analysis_mode = "standard_pando",
      stringsAsFactors = FALSE
    )
  )
  rownames(routing) <- NULL
  mode <- if (length(condition_types) && length(standard_types)) {
    "mixed_pando"
  } else if (length(condition_types)) {
    "condition_grn"
  } else {
    "standard_pando"
  }
  answer <- list(
    schema_version = "regcompass_celltype_routed_pando_v1",
    analysis_mode = mode,
    cell_type_analysis_mode = routing,
    condition_coefficients_calculated = length(condition_fits) > 0L,
    pando_fit_schema = if (length(condition_fits)) {
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
    } else NA_character_,
    pando_installed_version = as.character(utils::packageVersion("Pando")),
    pando_grn_data = if (is.null(condition_result)) NULL else
      condition_result$pando_grn_data,
    standard_pando_objects = standard_objects,
    condition_grn_fits = condition_fits,
    target_metabolic_genes = unique(unlist(
      lapply(results, `[[`, "target_metabolic_genes"), use.names = FALSE
    )),
    condition_fit_status = .rc_bind_pando_field(results, "condition_fit_status"),
    pando_network_index = if (is.null(condition_result)) data.frame() else
      condition_result$pando_network_index,
    pando_fit_diagnostics = if (is.null(condition_result)) data.frame() else
      condition_result$pando_fit_diagnostics,
    pando_execution_summary = if (is.null(condition_result)) list() else
      condition_result$pando_execution_summary,
    tf_peak_gene_universal = .rc_bind_pando_field(
      results, "tf_peak_gene_universal"
    ),
    tf_peak_gene_condition_all = .rc_bind_pando_field(
      results, "tf_peak_gene_condition_all"
    ),
    tf_peak_gene_condition = .rc_bind_pando_field(
      results, "tf_peak_gene_condition"
    ),
    tf_peak_gene_condition_effect_all = if (is.null(condition_result)) {
      data.frame()
    } else condition_result$tf_peak_gene_condition_effect_all,
    tf_peak_gene_condition_effect = if (is.null(condition_result)) {
      data.frame()
    } else condition_result$tf_peak_gene_condition_effect,
    paired_cell_ids = unique(as.character(paired_meta$cell_id)),
    paired_cell_metadata = paired_meta,
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      cell_type_routing =
        "condition GRN for at least two retained conditions; standard Pando otherwise",
      condition_effect_filter = "estimable and BH adjusted P below 0.05",
      standard_edge_filter =
        "adjusted P below 0.05 and absolute estimate above 0.01",
      projection =
        "paired-cell TF-by-ATAC before exact SuperCell aggregation"
    ),
    group_cols = c(condition_col, celltype_col)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(answer$condition_fit_status,
    file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_all.tsv.gz"))
  .rc_write_tsv_gz(answer$tf_peak_gene_condition,
    file.path(outdir, "pando_tf_peak_gene_active.tsv.gz"))
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_project_pando_by_celltype <- function(
    grn_result, membership, unit_meta, genes,
    condition_col, celltype_col, rna_assay, atac_assay) {
  template <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  projection <- reliability <- template
  coverage <- list()
  origins <- schemas <- projection_names <- policies <- character()
  if (length(grn_result$condition_grn_fits)) {
    part <- .rc_condition_pando_projection(
      grn_result, membership, unit_meta, genes
    )
    projection <- .rc_overlay_projection(
      projection, part$projection, "Condition-Pando"
    )
    reliability <- .rc_overlay_projection(
      reliability, part$reliability, "Condition-Pando reliability"
    )
    coverage[[length(coverage) + 1L]] <- part$coverage
    origins <- c(origins, part$origin)
    schemas <- c(schemas, part$pando_schema)
    projection_names <- c(projection_names, part$projection_name)
    policies <- c(policies, part$nonestimable_policy)
  }
  if (length(grn_result$standard_pando_objects)) {
    standard <- .rc_standard_pando_projection(
      grn_result, membership, unit_meta, condition_col, celltype_col,
      rna_assay = rna_assay, atac_assay = atac_assay,
      target_genes = genes
    )
    projection <- .rc_overlay_projection(
      projection, standard$projection, "Standard-Pando"
    )
    reliability <- .rc_overlay_projection(
      reliability, standard$reliability, "Standard-Pando reliability"
    )
    coverage[[length(coverage) + 1L]] <- standard$coverage
    origins <- c(origins, standard$projection_origin)
    schemas <- c(schemas, "standard_pando_network")
    projection_names <- c(projection_names, "standard_pando_full_fit")
    policies <- c(policies, "not_applicable_standard_pando")
  }
  if (!length(origins)) stop("No Pando projection route is available.", call. = FALSE)
  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = paste(unique(origins), collapse = ";"),
    pando_schema = paste(unique(schemas), collapse = ";"),
    projection_name = paste(unique(projection_names), collapse = ";"),
    nonestimable_policy = paste(unique(policies), collapse = ";"),
    condition_coefficients_calculated =
      length(grn_result$condition_grn_fits) > 0L,
    cell_type_analysis_mode = grn_result$cell_type_analysis_mode
  )
}
