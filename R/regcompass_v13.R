# Integrated two-stage RegCompass v1.3 workflow.

.rc_rbind_fill <- function(items) {
  items <- items[vapply(items, function(x) is.data.frame(x) && nrow(x), logical(1))]
  if (!length(items)) return(data.frame())
  columns <- unique(unlist(lapply(items, colnames), use.names = FALSE))
  items <- lapply(items, function(x) {
    for (column in setdiff(columns, colnames(x))) x[[column]] <- NA
    x[, columns, drop = FALSE]
  })
  out <- do.call(rbind, items); rownames(out) <- NULL; out
}

.rc_cbind_reaction_matrices <- function(items) {
  items <- items[vapply(items, function(x) !is.null(dim(x)) && ncol(x) > 0L, logical(1))]
  if (!length(items)) return(matrix(numeric(), 0L, 0L))
  reactions <- unique(unlist(lapply(items, rownames), use.names = FALSE))
  output <- lapply(items, function(x) {
    x <- as.matrix(x)
    aligned <- matrix(NA_real_, length(reactions), ncol(x), dimnames = list(reactions, colnames(x)))
    aligned[rownames(x), ] <- x
    aligned
  })
  answer <- do.call(cbind, output)
  if (anyDuplicated(colnames(answer))) stop("Duplicated metacell IDs across completed strata.", call. = FALSE)
  answer
}

.rc_cbind_same_rows <- function(items, label) {
  items <- items[vapply(items, function(x) !is.null(dim(x)) && ncol(x) > 0L, logical(1))]
  if (!length(items)) stop("Completed strata did not return ", label, ".", call. = FALSE)
  reference <- rownames(items[[1L]])
  if (is.null(reference) || any(!vapply(items, function(x) identical(rownames(x), reference), logical(1)))) {
    stop(label, " feature rows differ across completed strata.", call. = FALSE)
  }
  answer <- do.call(cbind, items)
  if (anyDuplicated(colnames(answer))) stop("Duplicated metacell IDs across completed strata.", call. = FALSE)
  answer
}

.rc_merge_completed_layer1 <- function(records, gpr_table, calibration_args = list()) {
  rna <- .rc_cbind_same_rows(lapply(records, function(x) x$layer1$rna_metacell_logcpm), "RNA metacell logCPM")
  detection <- .rc_cbind_same_rows(lapply(records, function(x) x$layer1$rna_metacell_detection), "RNA metacell detection")
  worker_raw <- .rc_cbind_reaction_matrices(lapply(records, function(x) x$layer1$C_raw %||% x$layer1$reaction_capacity_L1))
  confidence <- .rc_cbind_reaction_matrices(lapply(records, function(x) {
    capacity <- x$layer1$C_raw %||% x$layer1$reaction_capacity_L1
    rc_layer2_confidence_matrix(x$layer1$reaction_confidence, capacity)
  }))
  unit_meta <- .rc_rbind_fill(lapply(records, function(x) x$layer1$unit_meta))
  id_columns <- intersect(c("pool_id", "metacell_id", "unit_id"), colnames(unit_meta))
  if (!length(id_columns)) stop("Merged unit metadata lack a metacell identifier.", call. = FALSE)
  unit_meta <- unit_meta[match(colnames(rna), as.character(unit_meta[[id_columns[[1L]]]])), , drop = FALSE]
  if (anyNA(unit_meta[[id_columns[[1L]]]])) stop("Merged unit metadata are incomplete.", call. = FALSE)
  defaults <- list(
    gpr_table = gpr_table, unit_expression = rna, unit_detection = detection,
    unit_meta = unit_meta, stratum_col = NULL, gene_confidence = NULL,
    run_sensitivity = FALSE, bootstrap = FALSE, BPPARAM = FALSE
  )
  defaults[names(calibration_args)] <- NULL
  global <- do.call(rc_run_layer1_capacity, c(defaults, calibration_args))
  global$worker_C_raw <- worker_raw
  global$reaction_confidence <- confidence
  global$unit_meta <- unit_meta
  global$metacell_meta <- unit_meta
  global$rna_metacell_logcpm <- rna
  global$rna_metacell_detection <- detection
  global$global_calibration <- TRUE
  global$calibration_scope <- "all_completed_metacells"
  global$upstream_strata <- .rc_rbind_fill(lapply(records, function(x) x$status))
  global
}

.rc_merge_completed_modules <- function(records) {
  fields <- c(
    "sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  output <- lapply(fields, function(field) .rc_rbind_fill(lapply(records, function(x) x$grn_meta_modules[[field]])))
  names(output) <- fields
  output$schema_version <- "regcompass_global_union_meta_module_v1.3"
  output$crossref_maps <- records[[1L]]$grn_meta_modules$crossref_maps
  output$completed_strata <- .rc_rbind_fill(lapply(records, function(x) x$status))
  output
}

.rc_stop_upstream_stage <- function(BPPARAM) {
  .rc_stop_bpparam(BPPARAM)
  invisible(gc(verbose = FALSE))
}

#' Run the integrated strict-stratum RegCompass workflow
#'
#' Every retained condition x sample x cell-type stratum completes metacell,
#' LinkPeaks/Layer 1, Pando, and meta-module inference inside one upstream worker.
#' A hard barrier then merges all successful strata, recalibrates raw capacities
#' across all metacells, builds one global-union GEM per medium, and starts a
#' separate metacell-level LP parallel stage.
#' @export
rc_run_regcompass <- function(object, gem, outdir, pfm, genome,
                               fragment_files = NULL,
                               sample_col = "sample_id",
                               condition_col = "condition",
                               celltype_col = "cell_type",
                               rna_assay = "RNA",
                               atac_assay = "ATAC",
                               model_mode = c("meta_module_gem", "full_gem"),
                               medium_scenarios = NULL,
                               metacell_args = list(),
                               pando_args = list(),
                               calibration_args = list(bootstrap = FALSE),
                               layer2_args = list(),
                               BPPARAM_upstream = NULL,
                               BPPARAM_layer2 = NULL) {
  model_mode <- match.arg(model_mode)
  if (!inherits(object, "Seurat")) stop("`object` must inherit from Seurat.", call. = FALSE)
  if (is.null(gem$gpr_table)) stop("`gem` must contain `gpr_table`.", call. = FALSE)
  strict_cols <- .rc_strict_stratum_cols(sample_col, condition_col, celltype_col)
  missing <- setdiff(strict_cols, colnames(object@meta.data))
  if (length(missing)) stop("Missing strict-stratum metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  reserved_metacell <- intersect(names(metacell_args), c(
    "object", "gpr_table", "outdir", "fragment_files", "sample_col", "condition_col",
    "celltype_col", "rna_assay", "atac_assay", "BPPARAM_metacell", "BPPARAM_linkpeaks",
    "BPPARAM_layer1", "future_plan"
  ))
  if (length(reserved_metacell)) stop("`metacell_args` cannot override integrated fields: ", paste(reserved_metacell, collapse = ", "), call. = FALSE)
  reserved_pando <- intersect(names(pando_args), c(
    "metacell_object", "gem", "outdir", "pfm", "genome", "sample_col", "condition_col",
    "celltype_col", "single_cell_genes", "rna_assay", "atac_assay", "group_cols", "BPPARAM",
    "on_sample_error"
  ))
  if (length(reserved_pando)) stop("`pando_args` cannot override integrated fields: ", paste(reserved_pando, collapse = ", "), call. = FALSE)
  reserved_layer2 <- intersect(names(layer2_args), c(
    "layer1", "gem", "mode", "reaction_membership", "core_reactions", "medium_scenarios",
    "sample_col", "condition_col", "celltype_col", "unit", "BPPARAM"
  ))
  if (length(reserved_layer2)) stop("`layer2_args` cannot override integrated fields: ", paste(reserved_layer2, collapse = ", "), call. = FALSE)
  meta <- object@meta.data
  meta$.rc_stratum_id <- rc_make_stratum_id(meta, strict_cols)
  meta$.rc_cell_id <- rownames(meta)
  min_cells <- as.integer(metacell_args$min_cells_pre_metacell %||% 100L)
  counts <- table(meta$.rc_stratum_id)
  retained_ids <- names(counts[counts >= min_cells])
  if (!length(retained_ids)) stop("No strict stratum passes `min_cells_pre_metacell`.", call. = FALSE)
  cells_by_stratum <- split(meta$.rc_cell_id[meta$.rc_stratum_id %in% retained_ids], meta$.rc_stratum_id[meta$.rc_stratum_id %in% retained_ids])
  expected <- do.call(rbind, lapply(names(cells_by_stratum), function(id) {
    cells <- cells_by_stratum[[id]]
    values <- meta[match(cells[[1L]], meta$.rc_cell_id), strict_cols, drop = FALSE]
    data.frame(stratum_id = id, values, n_input_cells = length(cells), stringsAsFactors = FALSE, check.names = FALSE)
  }))
  .rc_write_tsv_gz(expected, file.path(outdir, "upstream_expected_strata.tsv.gz"))
  single_cell_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  run_one_stratum <- function(stratum_id) {
    cells <- cells_by_stratum[[stratum_id]]
    values <- meta[match(cells[[1L]], meta$.rc_cell_id), strict_cols, drop = FALSE]
    safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", stratum_id)
    stratum_dir <- file.path(outdir, "01_strata", safe_id)
    status <- data.frame(
      stratum_id = stratum_id, values, n_input_cells = length(cells),
      status = "pending", result_file = NA_character_, error_class = NA_character_,
      error_message = NA_character_, stringsAsFactors = FALSE, check.names = FALSE
    )
    answer <- tryCatch({
      subobject <- subset(object, cells = cells)
      metacell_defaults <- list(
        object = subobject, gpr_table = gem$gpr_table, outdir = stratum_dir,
        fragment_files = fragment_files, sample_col = sample_col, condition_col = condition_col,
        celltype_col = celltype_col, rna_assay = rna_assay, atac_assay = atac_assay,
        BPPARAM_metacell = FALSE, BPPARAM_linkpeaks = FALSE, BPPARAM_layer1 = FALSE,
        future_plan = "sequential"
      )
      metacell_defaults[names(metacell_args)] <- NULL
      layer1 <- do.call(rc_run_regcompass_multiome_metacell, c(metacell_defaults, metacell_args))
      retained <- as.character(layer1$metacell_meta$metacell_id)
      if (!length(retained)) stop("Stratum produced no retained metacells.", call. = FALSE)
      metacell_object <- rc_load_metacell_object_from_run(
        stratum_dir, retained_metacell_ids = retained,
        rna_assay = rna_assay, atac_assay = atac_assay
      )
      pando_defaults <- list(
        metacell_object = metacell_object, gem = gem,
        outdir = file.path(stratum_dir, "04_pando_meta_modules"), pfm = pfm, genome = genome,
        sample_col = sample_col, condition_col = condition_col, celltype_col = celltype_col,
        group_cols = strict_cols, single_cell_genes = single_cell_genes,
        rna_assay = rna_assay, atac_assay = atac_assay, BPPARAM = FALSE,
        on_sample_error = "stop"
      )
      pando_defaults[names(pando_args)] <- NULL
      modules <- do.call(rc_run_pando_meta_modules, c(pando_defaults, pando_args))
      result_file <- file.path(stratum_dir, "completed_stratum.rds")
      status$status <- "ok"; status$result_file <- result_file
      saveRDS(list(status = status, layer1 = layer1, grn_meta_modules = modules), result_file)
      status
    }, error = function(e) {
      status$status <- "failed"
      status$error_class <- class(e)[[1L]]
      status$error_message <- conditionMessage(e)
      status
    })
    answer
  }
  upstream_status <- tryCatch(
    rc_parallel_lapply(names(cells_by_stratum), run_one_stratum, BPPARAM = BPPARAM_upstream),
    finally = .rc_stop_upstream_stage(BPPARAM_upstream)
  )
  upstream_status <- do.call(rbind, upstream_status)
  rownames(upstream_status) <- NULL
  .rc_write_tsv_gz(upstream_status, file.path(outdir, "upstream_stratum_status.tsv.gz"))
  failed <- upstream_status$status != "ok" | is.na(upstream_status$result_file) | !file.exists(upstream_status$result_file)
  barrier_ids_match <- !anyDuplicated(upstream_status$stratum_id) &&
    setequal(as.character(upstream_status$stratum_id), as.character(expected$stratum_id))
  if (any(failed) || nrow(upstream_status) != nrow(expected) || !barrier_ids_match) {
    stop("Global barrier not reached: every retained stratum must finish Meta cell, LinkPeaks/Layer 1, Pando, and meta-module inference before downstream processing. See upstream_stratum_status.tsv.gz.", call. = FALSE)
  }
  records <- lapply(as.character(upstream_status$result_file), readRDS)
  layer1 <- .rc_merge_completed_layer1(records, gem$gpr_table, calibration_args)
  modules <- .rc_merge_completed_modules(records)
  if (!nrow(modules$core_gene_reaction) || !nrow(modules$reaction_membership)) {
    stop("Global barrier completed, but no meta-module reactions were available for union GEM construction.", call. = FALSE)
  }
  hard_core <- modules$core_gene_reaction
  if ("is_core" %in% colnames(hard_core)) hard_core <- hard_core[hard_core$is_core %in% TRUE, , drop = FALSE]
  default_targets <- unique(as.character(hard_core$reaction_id))
  saveRDS(layer1, file.path(outdir, "global_layer1_calibrated.rds"))
  saveRDS(modules, file.path(outdir, "global_meta_modules.rds"))
  layer2_defaults <- list(
    layer1 = layer1, gem = gem, target_reactions = default_targets,
    medium_scenarios = medium_scenarios, mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) modules$reaction_membership else NULL,
    core_reactions = if (identical(model_mode, "meta_module_gem")) modules$core_gene_reaction else NULL,
    unit = "metacell", sample_col = sample_col, condition_col = condition_col,
    celltype_col = celltype_col, BPPARAM = BPPARAM_layer2
  )
  layer2_defaults[names(layer2_args)] <- NULL
  microcompass <- do.call(rc_run_microcompass, c(layer2_defaults, layer2_args))
  result <- list(
    schema_version = "regcompass_v1.3_global_shared_gem",
    model_mode = model_mode, upstream_status = upstream_status,
    layer1 = layer1, grn_meta_modules = modules, microcompass = microcompass,
    barrier = list(expected_strata = nrow(expected), completed_strata = nrow(upstream_status), passed = TRUE),
    params = list(sample_col = sample_col, condition_col = condition_col,
                  celltype_col = celltype_col, shared_gem = TRUE,
                  layer2_unit = "metacell")
  )
  saveRDS(result, file.path(outdir, "regcompass_v1.3_result.rds"))
  result
}
