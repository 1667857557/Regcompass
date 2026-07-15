.rc_bind_frames_fill <- function(x) {
  x <- x[vapply(x, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
  if (!length(x)) return(data.frame())
  columns <- unique(unlist(lapply(x, colnames), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(columns, colnames(z))
    for (column in missing) z[[column]] <- NA
    z[, columns, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

.rc_cbind_matrix_union <- function(x, fill = NA_real_) {
  x <- x[vapply(x, function(z) !is.null(dim(z)) && ncol(z) > 0L, logical(1))]
  if (!length(x)) return(matrix(numeric(), 0L, 0L))
  rows <- unique(unlist(lapply(x, rownames), use.names = FALSE))
  columns <- unlist(lapply(x, colnames), use.names = FALSE)
  if (anyNA(rows) || any(!nzchar(rows)) || anyNA(columns) || any(!nzchar(columns))) {
    stop("Matrices require non-empty row and column names.", call. = FALSE)
  }
  if (anyDuplicated(columns)) {
    duplicated_ids <- unique(columns[duplicated(columns)])
    stop("Duplicated metacell IDs across strata: ", paste(utils::head(duplicated_ids, 10L), collapse = ", "), call. = FALSE)
  }
  out <- matrix(fill, nrow = length(rows), ncol = length(columns), dimnames = list(rows, columns))
  offset <- 0L
  for (matrix_in in x) {
    column_index <- seq.int(offset + 1L, offset + ncol(matrix_in))
    out[rownames(matrix_in), column_index] <- as.matrix(matrix_in)
    offset <- offset + ncol(matrix_in)
  }
  out
}

.rc_phase_bpparam <- function(workers = NULL, backend = c("auto", "serial", "snow", "multicore")) {
  backend <- match.arg(backend)
  if (identical(backend, "serial")) return(FALSE)
  param <- rc_default_bpparam(workers = workers, backend = backend)
  param %||% FALSE
}

.rc_release_bpparam <- function(param) {
  if (!identical(param, FALSE) && !is.null(param) && requireNamespace("BiocParallel", quietly = TRUE)) {
    try(BiocParallel::bpstop(param), silent = TRUE)
  }
  invisible(gc(verbose = FALSE))
}

.rc_run_regcompass_stratum <- function(object, group_id, group_cols, gem, outdir, pfm, genome,
                                        fragment_files = NULL, sample_col = "sample_id",
                                        condition_col = "condition", celltype_col = "cell_type",
                                        rna_assay = "RNA", atac_assay = "ATAC",
                                        metacell_args = list(), layer1_args = list(), pando_args = list()) {
  meta <- object@meta.data
  ids <- rc_make_stratum_id(meta, group_cols)
  cells <- rownames(meta)[ids == group_id]
  if (!length(cells)) stop("No cells found for stratum: ", group_id, call. = FALSE)
  one <- subset(object, cells = cells)
  stratum_dir <- file.path(outdir, gsub("[^A-Za-z0-9_.-]+", "_", group_id))
  dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)

  capacity_params <- list(
    promiscuity_mode = layer1_args$promiscuity_mode %||% "sqrt",
    and_method = layer1_args$and_method %||% "boltzmann",
    tau = layer1_args$tau %||% 0.20,
    or_method = "sum_sqrtK"
  )

  metacell_defaults <- list(
    object = one,
    outdir = file.path(stratum_dir, "01_metacells"),
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = TRUE,
    require_fragment_aggregation = TRUE,
    fragment_aggregation_backend = "regcompass",
    BPPARAM = FALSE,
    on_stratum_error = "stop"
  )
  reserved <- intersect(names(metacell_args), names(metacell_defaults))
  reserved <- setdiff(reserved, c("rna_reduction", "atac_reduction", "rna_dims", "atac_dims", "gamma",
                                  "seed", "min_cells_per_stratum", "min_metacell_size",
                                  "min_metacells_per_stratum", "adaptive_gamma", "label_col",
                                  "bgzip_path", "tabix_path", "fragment_nb_cl", "overwrite"))
  if (length(reserved)) stop("`metacell_args` cannot override workflow fields: ", paste(reserved, collapse = ", "), call. = FALSE)
  metacell_defaults[names(metacell_args)] <- NULL
  metacells <- do.call(rc_make_supercell2_metacells, c(metacell_defaults, metacell_args))

  minimum_metacells <- max(as.integer(c(
    layer1_args$min_metacells_for_linkpeaks %||% 10L,
    pando_args$min_metacells %||% 10L
  )), na.rm = TRUE)
  if (nrow(metacells$metacell_meta) < minimum_metacells) {
    stop("Stratum `", group_id, "` produced fewer than ", minimum_metacells, " metacells.", call. = FALSE)
  }
  metacell_object <- rc_load_or_merge_metacell_objects(
    metacells$metacell_objects,
    fragment_manifest = metacells$fragment_manifest,
    metacell_meta = metacells$metacell_meta,
    fragment_files = metacells$fragment_files,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = TRUE
  )

  linkpeaks_args <- layer1_args$linkpeaks_args %||% list()
  linkpeaks_args$genome <- linkpeaks_args$genome %||% genome
  linkpeaks_args$BPPARAM <- FALSE
  layer1_args$linkpeaks_args <- NULL
  layer1_args$BPPARAM <- NULL
  layer1_args$link_stratum_cols <- NULL
  layer1_args$stratum_col <- NULL
  layer1_args$min_metacells_for_linkpeaks <- NULL
  layer1_defaults <- list(
    gpr_table = gem$gpr_table,
    rna_metacell_counts = metacells$rna_counts,
    metacell_meta = metacells$metacell_meta,
    atac_metacell_counts = metacells$atac_counts,
    metacell_seurat = metacell_object,
    force_metacell_relink = TRUE,
    allow_supplied_links = FALSE,
    link_stratum_cols = group_cols,
    min_metacells_for_linkpeaks = minimum_metacells,
    metabolic_genes = gem$metabolic_genes %||% rc_metabolic_gpr_genes(gem$gpr_table),
    linkpeaks_args = linkpeaks_args,
    stratum_col = "stratum_id",
    BPPARAM = FALSE
  )
  layer1_defaults[names(layer1_args)] <- NULL
  layer1 <- do.call(rc_run_layer1_from_metacells, c(layer1_defaults, layer1_args))
  layer1$strict_group_id <- group_id

  pando_args$BPPARAM <- NULL
  pando_infer_args <- pando_args$pando_infer_args %||% list()
  pando_infer_args$parallel <- FALSE
  pando_args$pando_infer_args <- pando_infer_args
  pando_args$group_cols <- NULL
  pando_args$sample_col <- NULL
  pando_args$condition_col <- NULL
  pando_args$celltype_col <- NULL
  pando_defaults <- list(
    metacell_object = metacell_object,
    gem = gem,
    outdir = file.path(stratum_dir, "02_pando_meta_modules"),
    pfm = pfm,
    genome = genome,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    group_cols = group_cols,
    single_cell_genes = rownames(.rc_get_assay_counts(one, rna_assay)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_metacells = minimum_metacells,
    BPPARAM = FALSE,
    on_sample_error = "stop"
  )
  pando_defaults[names(pando_args)] <- NULL
  meta_modules <- do.call(rc_run_pando_meta_modules, c(pando_defaults, pando_args))

  artifact <- list(
    schema_version = "regcompass_stratum_v2",
    group_id = group_id,
    group_cols = group_cols,
    capacity_params = capacity_params,
    layer1 = layer1,
    grn_meta_modules = meta_modules,
    metacell_meta = metacells$metacell_meta,
    metacell_dir = stratum_dir
  )
  artifact_file <- file.path(stratum_dir, "stratum_result.rds")
  saveRDS(artifact, artifact_file)
  list(group_id = group_id, status = "ok", artifact_file = artifact_file,
       n_cells = length(cells), n_metacells = nrow(metacells$metacell_meta),
       error_class = NA_character_, error_message = NA_character_)
}

.rc_merge_stratum_layer1 <- function(artifacts, gem, single_cell_genes,
                                      sample_col, condition_col, celltype_col) {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  expression_list <- lapply(artifacts, function(x) x$layer1$rna_metacell_logcpm)
  if (any(vapply(expression_list, is.null, logical(1)))) {
    stop("Every upstream artifact must contain metacell RNA logCPM for global capacity recomputation.", call. = FALSE)
  }
  expression_list <- lapply(expression_list, function(x) {
    keep <- tolower(rownames(x)) %in% gpr_genes
    as.matrix(x[keep, , drop = FALSE])
  })
  common_genes <- Reduce(intersect, lapply(expression_list, rownames))
  common_genes <- rownames(expression_list[[1L]])[rownames(expression_list[[1L]]) %in% common_genes]
  if (!length(common_genes)) stop("No common GPR genes remain across upstream metacell artifacts.", call. = FALSE)
  expression_list <- lapply(expression_list, function(x) x[common_genes, , drop = FALSE])
  rna_logcpm <- .rc_cbind_matrix_union(expression_list)

  parameter_keys <- vapply(artifacts, function(x) {
    params <- x$capacity_params
    if (is.null(params)) return(NA_character_)
    paste(params$promiscuity_mode, params$and_method, params$tau, params$or_method, sep = "|")
  }, character(1))
  if (anyNA(parameter_keys) || length(unique(parameter_keys)) != 1L) {
    stop("Upstream artifacts must use one identical Layer 1 capacity parameter set.", call. = FALSE)
  }
  capacity_params <- artifacts[[1L]]$capacity_params
  global_gene_score <- rc_gene_score(rna_logcpm)
  C_raw <- rc_reaction_capacity(
    parsed,
    global_gene_score,
    promiscuity_mode = capacity_params$promiscuity_mode,
    tau = capacity_params$tau,
    and_method = capacity_params$and_method,
    or_method = capacity_params$or_method,
    BPPARAM = FALSE
  )
  calibrated <- rc_q95_calibrate(
    C_raw,
    bootstrap = FALSE,
    BPPARAM = FALSE,
    unit_meta = NULL,
    stratum_col = NULL
  )

  confidence_list <- lapply(artifacts, function(x) {
    rc_layer2_confidence_matrix(x$layer1$reaction_confidence, x$layer1$C_raw)
  })
  reaction_confidence <- .rc_cbind_matrix_union(confidence_list)
  reaction_confidence <- rc_align_layer2_evidence(
    reaction_confidence,
    rownames(C_raw),
    NA_real_
  )
  missing_confidence_units <- setdiff(colnames(C_raw), colnames(reaction_confidence))
  if (length(missing_confidence_units)) {
    stop("Reaction confidence is missing metacells after global alignment.", call. = FALSE)
  }
  reaction_confidence <- reaction_confidence[, colnames(C_raw), drop = FALSE]

  unit_meta <- .rc_bind_frames_fill(lapply(artifacts, function(x) x$layer1$unit_meta))
  id_col <- if ("pool_id" %in% colnames(unit_meta)) "pool_id" else if ("metacell_id" %in% colnames(unit_meta)) "metacell_id" else NULL
  if (is.null(id_col)) stop("Merged metacell metadata lacks pool_id/metacell_id.", call. = FALSE)
  unit_meta <- unit_meta[match(colnames(C_raw), as.character(unit_meta[[id_col]])), , drop = FALSE]
  if (anyNA(unit_meta[[id_col]])) stop("Merged metacell metadata are incomplete.", call. = FALSE)
  if (!"pool_id" %in% colnames(unit_meta)) unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  if (!"unit_id" %in% colnames(unit_meta)) unit_meta$unit_id <- unit_meta$pool_id
  unit_meta$stratum_id <- rc_make_stratum_id(unit_meta, c(condition_col, sample_col, celltype_col))

  list(
    schema_version = "regcompass_global_layer1_v2",
    C_or_raw = C_raw,
    C_raw = C_raw,
    reaction_capacity_L1 = C_raw,
    C_rel = calibrated$C_rel,
    reaction_confidence = reaction_confidence,
    q95_diagnostics = calibrated$Q,
    capacity_calibration_scope = "all_metacells_global_gene_score_and_reaction_q95",
    capacity_params = capacity_params,
    rna_metacell_logcpm = rna_logcpm,
    global_gene_score = global_gene_score,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, tolower(single_cell_genes)),
    parsed_gpr = parsed,
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "metacell",
    strict_stratum_cols = c(condition_col, sample_col, celltype_col)
  )
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  names_to_merge <- c("sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
                      "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
                      "reaction_membership", "meta_module_summary")
  out <- lapply(names_to_merge, function(name) {
    .rc_bind_frames_fill(lapply(artifacts, function(x) x$grn_meta_modules[[name]]))
  })
  names(out) <- names_to_merge
  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) core <- core[core$is_core %in% TRUE, , drop = FALSE]
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]
  membership_ids <- unique(as.character(out$reaction_membership$reaction_id))
  membership_ids <- membership_ids[!is.na(membership_ids) & nzchar(membership_ids)]
  if (!length(core_ids) || !length(membership_ids)) stop("No global meta-module core or membership reactions were produced.", call. = FALSE)
  out$global_core_reactions <- data.frame(
    sample_id = "global", module_id = "GLOBAL_UNION", reaction_id = core_ids,
    is_core = TRUE, stringsAsFactors = FALSE
  )
  out$global_reaction_membership <- data.frame(
    sample_id = "global", module_id = "GLOBAL_UNION", reaction_id = membership_ids,
    is_core = membership_ids %in% core_ids,
    inclusion_stage = ifelse(membership_ids %in% core_ids, "global_union_core", "global_union_member"),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_global_meta_module_v2"
  out$source_group_ids <- unique(unlist(lapply(artifacts, `[[`, "group_id"), use.names = FALSE))
  out
}
