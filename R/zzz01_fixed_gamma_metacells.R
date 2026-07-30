# Canonical fixed-gamma metacell construction.

.rc_make_condition_celltype_metacells_adaptive <-
  .rc_make_condition_celltype_metacells
.rc_make_condition_celltype_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list()) {
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  gamma <- metacell_args$gamma %||% 30L
  gamma <- as.integer(gamma)
  if (length(gamma) != 1L || is.na(gamma) || gamma < 1L) {
    stop("`metacell_args$gamma` must be one positive integer.", call. = FALSE)
  }
  metacell_args$gamma <- gamma
  metacell_args$depth_balance <- FALSE
  out <- .rc_make_condition_celltype_metacells_adaptive(
    object = object,
    outdir = outdir,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args
  )
  out$input_design$gamma <- gamma
  out$input_design$depth_balance <- FALSE
  out$input_design$depth_balance_policy <-
    "fixed_gamma_no_depth_restriction"
  out$input_design$extreme_depth_policy <-
    "diagnostic_only_no_top1_rejection"
  out
}

.rc_regcompass_step_metacells_adaptive <- rc_regcompass_step_metacells
rc_regcompass_step_metacells <- function(
    object, outdir,
    condition_col = "condition",
    celltype_col = "cell_type",
    cell_type = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    fragment_files = FALSE,
    metacell_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  if (!is.list(metacell_args)) {
    stop("`metacell_args` must be a list.", call. = FALSE)
  }
  metacell_args$gamma <- as.integer(metacell_args$gamma %||% 30L)
  metacell_args$depth_balance <- FALSE
  .rc_regcompass_step_metacells_adaptive(
    object = object,
    outdir = outdir,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    metacell_args = metacell_args,
    progress = progress
  )
}

.rc_make_supercell2_metacells_adaptive <- rc_make_supercell2_metacells
rc_make_supercell2_metacells <- function(
    object,
    outdir,
    sample_col = NULL,
    stratum_col = NULL,
    condition_col = "condition",
    celltype_col = "cell_type",
    state_col = NULL,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    rna_reduction = "pca",
    atac_reduction = "lsi",
    rna_dims = 1:30,
    atac_dims = 2:30,
    gamma = 30L,
    seed = 12345L,
    min_cells_per_stratum = 100L,
    min_metacell_size = 20L,
    min_metacells_per_stratum = 2L,
    depth_balance = FALSE,
    label_col = NULL,
    fragment_files = FALSE,
    bgzip_path = "bgzip",
    tabix_path = "tabix",
    fragment_nb_cl = 1L,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = FALSE,
    require_fragment_aggregation = FALSE,
    fragment_aggregation_backend = c("none", "regcompass", "supercell"),
    overwrite = FALSE,
    BPPARAM = NULL,
    on_stratum_error = c("record", "stop"),
    call_peaks_from_fragments = FALSE,
    macs2_path = NULL,
    peak_calling_effective_genome_size = NULL,
    peak_calling_args = list()) {
  on_stratum_error <- match.arg(on_stratum_error)
  fragment_aggregation_backend <- match.arg(fragment_aggregation_backend)
  if (is.null(stratum_col)) stratum_col <- sample_col
  if (!is.character(stratum_col) || length(stratum_col) != 1L ||
      is.na(stratum_col) || !nzchar(trimws(stratum_col))) {
    stop("`stratum_col` must identify one metadata column.", call. = FALSE)
  }
  controls <- c(
    gamma = gamma,
    min_cells_per_stratum = min_cells_per_stratum,
    min_metacell_size = min_metacell_size,
    min_metacells_per_stratum = min_metacells_per_stratum
  )
  if (any(!is.finite(controls)) || any(controls < 1) ||
      any(abs(controls - round(controls)) > sqrt(.Machine$double.eps))) {
    stop("Metacell controls must be positive integers.", call. = FALSE)
  }
  gamma <- as.integer(gamma)
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!isTRUE(save_counts)) {
    stop("`save_counts` must be TRUE.", call. = FALSE)
  }
  fragment_requested <- !identical(fragment_files, FALSE) &&
    !is.null(fragment_files)
  if (fragment_requested || isTRUE(save_fragments) ||
      isTRUE(require_fragment_aggregation) ||
      !identical(fragment_aggregation_backend, "none") ||
      isTRUE(call_peaks_from_fragments)) {
    stop(
      paste(
        "The fixed-gamma canonical path aggregates the existing RNA and ATAC",
        "count assays and requires fragment_files = FALSE, save_fragments =",
        "FALSE, require_fragment_aggregation = FALSE,",
        "fragment_aggregation_backend = 'none', and",
        "call_peaks_from_fragments = FALSE."
      ),
      call. = FALSE
    )
  }
  if (isTRUE(depth_balance)) {
    warning(
      "Ignoring `depth_balance = TRUE`; every stratum uses the same fixed gamma.",
      call. = FALSE
    )
  }
  .rc_require_supercell2()
  .rc_validate_supercell2_inputs(
    object,
    assays = c(rna_assay, atac_assay),
    reductions = c(rna_reduction, atac_reduction)
  )
  meta <- object@meta.data
  required <- unique(c(
    stratum_col, condition_col, celltype_col, state_col, label_col
  ))
  required <- required[
    !is.null(required) & !is.na(required) & nzchar(required)
  ]
  missing <- setdiff(required, colnames(meta))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  cell_order <- colnames(object)
  meta <- meta[match(cell_order, rownames(meta)), , drop = FALSE]
  if (anyNA(rownames(meta))) {
    stop("Seurat metadata cannot be aligned to cells.", call. = FALSE)
  }
  rna_depth <- Matrix::colSums(
    .rc_get_assay_counts_safe(object, rna_assay)[, cell_order, drop = FALSE]
  )
  atac_depth <- Matrix::colSums(
    .rc_get_assay_counts_safe(object, atac_assay)[, cell_order, drop = FALSE]
  )
  meta$.rc_rna_depth <- as.numeric(rna_depth[cell_order])
  meta$.rc_atac_depth <- as.numeric(atac_depth[cell_order])
  meta$cell_id <- cell_order
  strata <- as.character(meta[[stratum_col]])
  if (anyNA(strata) || any(!nzchar(trimws(strata)))) {
    stop("`stratum_col` contains missing or empty labels.", call. = FALSE)
  }
  groups <- split(meta$cell_id, strata)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  build_one <- function(key) {
    cells <- groups[[key]]
    one_meta <- meta[match(cells, meta$cell_id), , drop = FALSE]
    stratum_dir <- file.path(
      outdir, paste0("stratum=", .rc_safe_path_component(key))
    )
    required_files <- c(
      file.path(stratum_dir, "metacell_metadata.tsv.gz"),
      file.path(stratum_dir, "membership.tsv.gz"),
      file.path(stratum_dir, "rna_counts.rds"),
      file.path(stratum_dir, "atac_counts.rds")
    )
    if (isTRUE(save_metacell_object)) {
      required_files <- c(
        required_files, file.path(stratum_dir, "metacell_object.rds")
      )
    }
    if (!isTRUE(overwrite) && all(file.exists(required_files))) {
      meta_existing <- utils::read.delim(
        gzfile(file.path(stratum_dir, "metacell_metadata.tsv.gz")),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      return(list(
        status = data.frame(
          stratum_id = key,
          n_input_cells = length(cells),
          gamma = gamma,
          actual_metacells = nrow(meta_existing),
          status = "cached",
          output_dir = stratum_dir,
          stringsAsFactors = FALSE
        ),
        output_dir = stratum_dir
      ))
    }
    dir.create(file.path(stratum_dir, "qc"), recursive = TRUE,
               showWarnings = FALSE)
    if (length(cells) < as.integer(min_cells_per_stratum)) {
      status <- data.frame(
        stratum_id = key,
        n_input_cells = length(cells),
        gamma = gamma,
        actual_metacells = 0L,
        status = "skipped_below_min_cells",
        output_dir = NA_character_,
        stringsAsFactors = FALSE
      )
      return(list(status = status, output_dir = NA_character_))
    }
    prefix <- paste0(.rc_safe_path_component(key), "_MC_")
    args <- list(
      seurat = subset(object, cells = cells),
      assay = c(rna_assay, atac_assay),
      reduction = list(rna_reduction, atac_reduction),
      dims = list(as.integer(rna_dims), as.integer(atac_dims)),
      gamma = gamma,
      return.seurat = TRUE,
      prefixMC = prefix,
      seed = as.integer(seed) + match(key, names(groups)) - 1L
    )
    if (!is.null(label_col)) args$label <- label_col
    mc <- .rc_supercell2_scimplify_for_seurat(args)
    SuperCell::validate_metacell_output(
      mc,
      rna_assay = rna_assay,
      atac_assay = atac_assay,
      require_fragments = FALSE,
      require_complete_fragment_coverage = FALSE
    )
    rna_counts <- .rc_as_sparse(
      .rc_get_assay_counts_safe(mc, rna_assay)
    )
    atac_counts <- .rc_as_sparse(
      .rc_get_assay_counts_safe(mc, atac_assay)
    )
    mc_ids <- as.character(colnames(rna_counts))
    if (!identical(mc_ids, as.character(colnames(atac_counts)))) {
      stop("RNA and ATAC metacell IDs differ within stratum `", key, "`.",
           call. = FALSE)
    }
    if (length(mc_ids) < as.integer(min_metacells_per_stratum)) {
      status <- data.frame(
        stratum_id = key,
        n_input_cells = length(cells),
        gamma = gamma,
        actual_metacells = length(mc_ids),
        status = "skipped_below_min_metacells",
        output_dir = NA_character_,
        stringsAsFactors = FALSE
      )
      return(list(status = status, output_dir = NA_character_))
    }
    membership <- .rc_extract_supercell_membership(mc, cells, mc_ids)
    membership_index <- match(membership$cell_id, rownames(one_meta))
    if (anyNA(membership_index)) {
      membership_index <- match(membership$cell_id, one_meta$cell_id)
    }
    metadata_columns <- unique(c(
      stratum_col, condition_col, celltype_col, state_col, label_col
    ))
    metadata_columns <- metadata_columns[
      !is.null(metadata_columns) & !is.na(metadata_columns) &
        nzchar(metadata_columns)
    ]
    source_index <- match(membership$cell_id, one_meta$cell_id)
    if (anyNA(source_index)) {
      stop("Metacell membership cannot be aligned to source metadata.",
           call. = FALSE)
    }
    for (column in metadata_columns) {
      membership[[column]] <- as.character(one_meta[[column]][source_index])
    }
    mc_meta <- rc_build_metacell_metadata(membership)
    mc_meta$stratum_id <- key
    mc_meta$metacell_display_id <- as.character(mc_meta$metacell_id)
    mc_meta$effective_gamma <- gamma
    mc_meta$requested_gamma <- gamma
    mc_meta$fixed_gamma <- TRUE
    mc_meta$depth_balance <- FALSE
    mc_meta$depth_restriction_applied <- FALSE
    mc_meta$low_power_metacell <-
      !is.na(mc_meta$n_cells) & mc_meta$n_cells < min_metacell_size

    input_index <- match(membership$cell_id, one_meta$cell_id)
    rna_cutoff <- unname(stats::quantile(
      one_meta$.rc_rna_depth, 0.99, na.rm = TRUE
    ))
    atac_cutoff <- unname(stats::quantile(
      one_meta$.rc_atac_depth, 0.99, na.rm = TRUE
    ))
    extreme <- one_meta$.rc_rna_depth[input_index] > rna_cutoff |
      one_meta$.rc_atac_depth[input_index] > atac_cutoff
    extreme_count <- tapply(
      extreme, as.character(membership$metacell_id), sum
    )
    mc_meta$extreme_depth_cells <- as.integer(
      extreme_count[match(mc_meta$metacell_id, names(extreme_count))]
    )
    mc_meta$extreme_depth_policy <-
      "diagnostic_only_no_top1_rejection"
    mc_meta$rna_total_umi <- as.numeric(
      Matrix::colSums(rna_counts)[mc_meta$metacell_id]
    )
    mc_meta$atac_total_fragments <- as.numeric(
      Matrix::colSums(atac_counts)[mc_meta$metacell_id]
    )
    mc_meta$depth_balance_policy <-
      "fixed_gamma_no_depth_restriction"

    if (isTRUE(save_metacell_object)) {
      saveRDS(mc, file.path(stratum_dir, "metacell_object.rds"))
    }
    saveRDS(rna_counts, file.path(stratum_dir, "rna_counts.rds"))
    saveRDS(atac_counts, file.path(stratum_dir, "atac_counts.rds"))
    .rc_write_tsv_gz(
      membership, file.path(stratum_dir, "membership.tsv.gz")
    )
    .rc_write_tsv_gz(
      mc_meta, file.path(stratum_dir, "metacell_metadata.tsv.gz")
    )
    qc <- data.frame(
      stratum_id = key,
      n_input_cells = length(cells),
      n_metacells = length(mc_ids),
      gamma = gamma,
      fixed_gamma = TRUE,
      depth_balance = FALSE,
      top1_depth_rejection = FALSE,
      max_top1_depth_cells_per_metacell = max(
        mc_meta$extreme_depth_cells, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
    .rc_write_tsv_gz(qc, file.path(stratum_dir, "qc", "metacell_qc.tsv.gz"))
    list(
      status = data.frame(
        stratum_id = key,
        n_input_cells = length(cells),
        gamma = gamma,
        actual_metacells = length(mc_ids),
        status = "ok",
        output_dir = stratum_dir,
        stringsAsFactors = FALSE
      ),
      output_dir = stratum_dir
    )
  }

  build_safe <- function(key) {
    tryCatch(
      build_one(key),
      error = function(error) {
        if (identical(on_stratum_error, "stop")) stop(error)
        list(
          status = data.frame(
            stratum_id = key,
            n_input_cells = length(groups[[key]]),
            gamma = gamma,
            actual_metacells = 0L,
            status = "failed",
            output_dir = NA_character_,
            error_message = conditionMessage(error),
            stringsAsFactors = FALSE
          ),
          output_dir = NA_character_
        )
      }
    )
  }
  results <- rc_parallel_lapply(names(groups), build_safe, BPPARAM = BPPARAM)
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  .rc_write_tsv_gz(
    status, file.path(outdir, "metacell_stratum_status.tsv.gz")
  )
  directories <- vapply(results, `[[`, character(1), "output_dir")
  directories <- directories[!is.na(directories) & nzchar(directories)]
  if (!length(directories)) {
    stop("All fixed-gamma metacell strata failed or were skipped.",
         call. = FALSE)
  }
  out <- rc_import_supercell2_metacells(
    directories,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    sample_col = stratum_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    require_fragments = FALSE
  )
  out$stratum_status <- status
  out$atac_count_source <- "aggregated_object_peak_counts"
  out$atac_peak_source <- "existing_object_peak_ranges"
  out$construction_policy <- list(
    gamma = gamma,
    gamma_scope = "global_fixed_across_all_strata",
    depth_balance = FALSE,
    top1_depth_rejection = FALSE,
    extreme_depth_policy = "diagnostic_only_no_restriction"
  )
  out
}
