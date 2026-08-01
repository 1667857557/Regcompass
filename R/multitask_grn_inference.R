.rc_run_celltype_multitask_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_design_args = list(),
    multitask_args = list(),
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_celltype_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  on_celltype_error <- match.arg(on_celltype_error)
  species <- .rc_infer_gem_species(gem, species)
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 3 ||
      abs(min_cells - round(min_cells)) > sqrt(.Machine$double.eps)) {
    stop("`min_cells` must be one integer of at least 3.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)
  if (!is.logical(save_pando_objects) || length(save_pando_objects) != 1L ||
      is.na(save_pando_objects)) {
    stop("`save_pando_objects` must be TRUE or FALSE.", call. = FALSE)
  }
  for (value in list(pando_initiate_args, pando_motif_args, pando_design_args)) {
    if (!is.list(value)) {
      stop("Pando argument bundles must be lists.", call. = FALSE)
    }
  }
  .rc_assert_pando_nested_boundary(
    pando_motif_args,
    c("object", "pfm", "genome", "store_motif_positions"),
    "pando_motif_args"
  )
  multitask_args <- .rc_validate_multitask_grn_args(multitask_args)
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Install 1667857557/Pando_regcompass before GRN inference.",
         call. = FALSE)
  }
  pando_version <- tryCatch(
    utils::packageVersion("Pando"), error = function(error) NULL
  )
  if (is.null(pando_version) || pando_version < "1.1.3" ||
      !exists("prepare_grn_design", envir = asNamespace("Pando"),
              inherits = FALSE) ||
      !exists("validate_grn_design", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "Multitask GRN inference requires Pando >= 1.1.3 from ",
      "1667857557/Pando_regcompass.", call. = FALSE
    )
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for multitask GRN inference.",
         call. = FALSE)
  }
  pando_install <- .rc_validate_pando_repository()
  motif_policy <- "user_supplied"
  if (is.null(pfm)) {
    pfm <- .rc_default_pando_motifs()
    motif_policy <- "Pando::motifs"
  }
  if (!length(pfm)) stop("`pfm` must be non-empty.", call. = FALSE)

  required_meta <- c(condition_col, celltype_col)
  missing <- setdiff(required_meta, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invalid_meta <- vapply(
    object@meta.data[, required_meta, drop = FALSE],
    function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))),
    logical(1)
  )
  if (any(invalid_meta)) {
    stop("Condition and cell-type metadata must be complete.", call. = FALSE)
  }
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop("Multitask Pando requires cell-type-shared ATAC TF-IDF.",
         call. = FALSE)
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"), recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_designs"), recursive = TRUE,
             showWarnings = FALSE)

  region_policy <- "user_supplied"
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
    region_policy <- if (identical(species, "human")) {
      paste(
        "union(Pando::phastConsElements20Mammals.UCSC.hg38,",
        "Pando::SCREEN.ccRE.UCSC.hg38)"
      )
    } else {
      "Pando::phastConsElements20Mammals.UCSC.hg38"
    }
  }

  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(
    toupper(.rc_mm_trim_unique(rna_genes)),
    toupper(.rc_mm_trim_unique(metabolic_genes))
  )
  target_genes <- .rc_mm_trim_unique(
    rna_genes[toupper(rna_genes) %in% target_upper]
  )
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }

  meta <- object@meta.data
  meta[[condition_col]] <- trimws(as.character(meta[[condition_col]]))
  meta[[celltype_col]] <- trimws(as.character(meta[[celltype_col]]))
  celltypes <- sort(unique(meta[[celltype_col]]))
  conditions <- sort(unique(meta[[condition_col]]))
  expected_groups <- expand.grid(
    condition_value = conditions,
    celltype_value = celltypes,
    stringsAsFactors = FALSE
  )
  names(expected_groups) <- c(condition_col, celltype_col)
  expected_groups$group_id <- rc_make_stratum_id(
    expected_groups, c(condition_col, celltype_col)
  )

  run_one_celltype <- function(celltype_value) {
    cells <- rownames(meta)[meta[[celltype_col]] == celltype_value]
    cell_meta <- meta[cells, , drop = FALSE]
    count <- table(factor(cell_meta[[condition_col]], levels = conditions))
    status <- data.frame(
      celltype_value = celltype_value,
      n_cells = length(cells),
      n_conditions = sum(count > 0),
      min_cells_per_condition = if (length(count)) min(count) else 0,
      n_structural_candidates = 0L,
      n_model_edges = 0L,
      n_active_condition_edges = 0L,
      status = "pending",
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    names(status)[1L] <- celltype_col
    if (any(count < min_cells)) {
      status$status <- "failed_too_few_cells"
      status$error_message <- paste0(
        "Every condition requires at least ", min_cells,
        " cells in cell type ", celltype_value, "."
      )
      return(list(
        status = status, candidates = data.frame(), global = data.frame(),
        condition = data.frame(), diagnostics = data.frame(), design = NULL
      ))
    }

    one <- tryCatch({
      obj <- subset(object, cells = cells)
      filtered <- .rc_drop_zero_count_atac_features(
        obj, atac_assay, paste0("Multitask Pando cell type ", celltype_value)
      )
      obj <- filtered$object
      init_defaults <- list(
        object = obj, peak_assay = atac_assay, rna_assay = rna_assay
      )
      init_defaults[names(pando_initiate_args)] <- NULL
      grn <- do.call(
        Pando::initiate_grn, c(init_defaults, pando_initiate_args)
      )
      motif_defaults <- list(
        object = grn,
        pfm = pfm,
        genome = genome,
        store_motif_positions = FALSE
      )
      motif_defaults[names(pando_motif_args)] <- NULL
      grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))
      design_defaults <- list(object = grn, genes = target_genes)
      design_defaults[names(pando_design_args)] <- NULL
      design <- do.call(
        Pando::prepare_grn_design,
        c(design_defaults, pando_design_args)
      )
      Pando::validate_grn_design(design)
      if (!identical(design$schema_version, "pando_grn_design_v2")) {
        stop("Pando did not return the required version-2 GRN design contract.",
             call. = FALSE)
      }
      rna <- .rc_pando_assay_data(obj, rna_assay)
      atac <- .rc_pando_assay_data(obj, atac_assay)
      fit <- .rc_fit_multitask_celltype_grn(
        design = design,
        rna = rna,
        atac = atac,
        meta = cell_meta,
        condition_col = condition_col,
        multitask_args = multitask_args
      )
      candidates <- design$candidate_edges
      candidates[[celltype_col]] <- celltype_value
      candidates$edge_universe_id <- design$design_fingerprint
      candidates <- candidates[, c(
        celltype_col, "edge_universe_id",
        setdiff(colnames(candidates), c(celltype_col, "edge_universe_id"))
      ), drop = FALSE]
      add_celltype <- function(value) {
        if (!nrow(value)) return(value)
        value[[celltype_col]] <- celltype_value
        value$edge_universe_id <- design$design_fingerprint
        value[, c(
          celltype_col, "edge_universe_id",
          setdiff(colnames(value), c(celltype_col, "edge_universe_id"))
        ), drop = FALSE]
      }
      fit$global <- add_celltype(fit$global)
      fit$condition <- add_celltype(fit$condition)
      fit$diagnostics <- add_celltype(fit$diagnostics)
      if (isTRUE(save_pando_objects)) {
        safe <- gsub("[^A-Za-z0-9_.-]+", "_", celltype_value)
        saveRDS(grn, file.path(outdir, "pando_objects", paste0(safe, ".rds")))
        saveRDS(
          design,
          file.path(outdir, "pando_designs", paste0(safe, ".rds"))
        )
      }
      list(
        candidates = candidates,
        global = fit$global,
        condition = fit$condition,
        diagnostics = fit$diagnostics,
        design = design,
        peak_diagnostics = filtered$diagnostics
      )
    }, error = function(error) error)

    if (inherits(one, "error")) {
      status$status <- "failed"
      status$error_class <- class(one)[[1L]]
      status$error_message <- conditionMessage(one)
      if (identical(on_celltype_error, "stop")) stop(one)
      return(list(
        status = status, candidates = data.frame(), global = data.frame(),
        condition = data.frame(), diagnostics = data.frame(), design = NULL
      ))
    }
    status$status <- "ok"
    status$n_structural_candidates <- nrow(one$candidates)
    status$n_model_edges <- nrow(one$global)
    status$n_active_condition_edges <- sum(
      one$condition$active_edge %in% TRUE, na.rm = TRUE
    )
    c(list(status = status), one)
  }

  results <- rc_parallel_lapply(
    celltypes, run_one_celltype, BPPARAM = BPPARAM
  )
  celltype_status <- do.call(rbind, lapply(results, `[[`, "status"))
  failed <- celltype_status$status != "ok"
  if (any(failed)) {
    stop(
      "Every cell-type multitask GRN must complete successfully. Failed: ",
      paste(celltype_status[[celltype_col]][failed], collapse = "; "),
      call. = FALSE
    )
  }
  bind <- function(name) {
    .rc_bind_frames_fill(lapply(results, `[[`, name))
  }
  candidates <- bind("candidates")
  global <- bind("global")
  condition_all <- bind("condition")
  diagnostics <- bind("diagnostics")
  if (!nrow(condition_all)) {
    stop("Multitask GRN inference produced no fitted condition edges.",
         call. = FALSE)
  }
  condition_all$condition <- as.character(condition_all$condition)
  names(condition_all)[names(condition_all) == "condition"] <- condition_col
  condition_all$group_id <- rc_make_stratum_id(
    condition_all, c(condition_col, celltype_col)
  )
  condition_all <- condition_all[, c(
    "group_id", condition_col, celltype_col,
    setdiff(colnames(condition_all), c("group_id", condition_col, celltype_col))
  ), drop = FALSE]
  significant <- condition_all[
    condition_all$active_edge %in% TRUE, , drop = FALSE
  ]
  if (!nrow(significant)) {
    stop("No bootstrap-stable multitask TF-peak-gene edges were available.",
         call. = FALSE)
  }

  condition_targets <- unique(significant[, c(
    "group_id", condition_col, celltype_col, "target"
  ), drop = FALSE])
  active_key <- paste(significant$group_id, significant$target, sep = "\001")
  condition_targets$n_active_edges <- as.integer(table(active_key)[paste(
    condition_targets$group_id, condition_targets$target, sep = "\001"
  )])
  condition_targets$evidence_definition <-
    "bootstrap_stability_selected_multitask_tf_peak_gene_target"

  status_rows <- lapply(seq_len(nrow(expected_groups)), function(i) {
    one <- expected_groups[i, , drop = FALSE]
    selected_cells <- meta[[condition_col]] == one[[condition_col]] &
      meta[[celltype_col]] == one[[celltype_col]]
    selected_edges <- significant[[condition_col]] == one[[condition_col]] &
      significant[[celltype_col]] == one[[celltype_col]]
    selected_targets <- condition_targets[[condition_col]] == one[[condition_col]] &
      condition_targets[[celltype_col]] == one[[celltype_col]]
    data.frame(
      one,
      n_cells = sum(selected_cells),
      n_active_edges = sum(selected_edges),
      n_active_targets = sum(selected_targets),
      status = "ok",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  group_status <- do.call(rbind, status_rows)
  rownames(group_status) <- NULL

  stability_columns <- intersect(c(
    "group_id", condition_col, celltype_col, "edge_universe_id",
    "edge_id", "tf", "region", "atac_feature_id", "target",
    "global_estimate", "condition_deviation", "effective_estimate",
    "selection_frequency", "sign_stability", "stability_weight",
    "stable_estimate", "active_edge", "sign_flip_flag", "cv_rsq",
    "bootstrap_method", "n_bootstrap_requested", "n_bootstrap_success"
  ), colnames(condition_all))
  stability <- condition_all[, stability_columns, drop = FALSE]

  .rc_mm_write_tsv_gz(
    candidates, file.path(outdir, "pando_tf_peak_gene_candidates.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    global, file.path(outdir, "pando_tf_peak_gene_global.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    condition_all, file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    significant, file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    condition_targets, file.path(outdir, "condition_target_genes.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    diagnostics, file.path(outdir, "target_model_diagnostics.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    stability, file.path(outdir, "bootstrap_stability_diagnostics.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    celltype_status, file.path(outdir, "pando_celltype_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    group_status, file.path(outdir, "pando_group_status.tsv.gz")
  )

  list(
    schema_version = "regcompass_multitask_grn_v2",
    grn_mode = "multitask_shared_backbone",
    target_metabolic_genes = target_genes,
    tf_peak_gene_candidates = candidates,
    tf_peak_gene_global = global,
    tf_peak_gene_condition_all = condition_all,
    tf_peak_gene_all = condition_all,
    tf_peak_gene_significant = significant,
    condition_target_genes = condition_targets,
    target_model_diagnostics = diagnostics,
    stability_diagnostics = stability,
    celltype_fit_status = celltype_status,
    group_status = group_status,
    group_cols = c(condition_col, celltype_col),
    pando_objects = if (isTRUE(save_pando_objects)) {
      file.path(outdir, "pando_objects")
    } else {
      NULL
    },
    pando_designs = if (isTRUE(save_pando_objects)) {
      file.path(outdir, "pando_designs")
    } else {
      NULL
    },
    normalization_policy = list(
      rna = "log-normalized RNA",
      atac = normalization,
      candidate_background =
        "one validated Pando structural design per cell type across conditions",
      residualization = "condition-centred target and TF-by-ATAC predictors",
      condition_balance = "equal total regression loss weight per condition",
      condition_parameterization =
        "global coefficient plus symmetric zero-sum condition deviations",
      stability = paste(
        "full-size nonparametric bootstrap with replacement within each",
        "condition; bootstrap samples are re-centred within condition and",
        "fitted at the full-data selected lambda"
      ),
      downstream_projection = paste(
        "bootstrap-stable TF-peak-target coefficients are projected with",
        "metacell ATAC only; condition-specific TF RNA is not reused"
      ),
      pando_motifs = motif_policy,
      pando_regions = region_policy,
      pando_design_schema = "pando_grn_design_v2"
    ),
    multitask_args = multitask_args,
    pando_installation = pando_install,
    pando_version = as.character(pando_version),
    species = species
  )
}
