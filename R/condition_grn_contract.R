# Authoritative Pando ConditionGRNFit integration.

.rc_reference_contrast <- function(beta, reference_condition) {
  beta <- as.matrix(beta)
  if (!is.character(reference_condition) ||
      length(reference_condition) != 1L ||
      !reference_condition %in% colnames(beta)) {
    stop("Reference condition is absent from the coefficient matrix.",
         call. = FALSE)
  }
  sweep(beta, 1L, beta[, reference_condition], "-")
}

.rc_extract_condition_grn_contract <- function(
    grn_object, condition_col, celltype_col,
    min_abs_estimate = 0, min_model_rsq = 0.1) {
  fits <- Pando::condition_grn_fit(
    grn_object, network_name = "regcompass_condition_grn"
  )
  if (inherits(fits, "ConditionGRNFit")) {
    fits <- list(fits)
  }
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), "ConditionGRNFit"))) {
    stop("Pando did not return complete ConditionGRNFit contracts.",
         call. = FALSE)
  }
  active_tol <- max(
    suppressWarnings(as.numeric(min_abs_estimate)), 1e-8, na.rm = TRUE
  )
  condition_rows <- list()
  universal_rows <- list()
  row_index <- 1L
  for (fit in fits) {
    if (!identical(fit$schema_version, "pando_condition_grn_fit_v2") ||
        !identical(
          fit$fit_engine, "shared_design_independent_elastic_net"
        ) ||
        !identical(
          fit$coefficient_scale,
          "pooled_cell_type_edge_and_target_standardized"
        )) {
      stop(
        "RegCompass requires a pooled-scale shared_design_independent ",
        "ConditionGRNFit v2.", call. = FALSE
      )
    }
    edge <- as.data.frame(fit$edge_table, stringsAsFactors = FALSE)
    transform <- as.data.frame(
      fit$predictor_transform, stringsAsFactors = FALSE
    )
    response <- as.data.frame(
      fit$response_transform, stringsAsFactors = FALSE
    )
    required_edge <- c("edge_id", "tf", "target", "region", "term")
    if (!all(required_edge %in% colnames(edge)) ||
        anyDuplicated(edge$edge_id)) {
      stop("ConditionGRNFit edge dictionary is invalid.", call. = FALSE)
    }
    if (!all(c("edge_id", "center", "scale") %in% colnames(transform)) ||
        anyDuplicated(transform$edge_id)) {
      stop("ConditionGRNFit predictor transform is invalid.",
           call. = FALSE)
    }
    transform <- transform[
      match(edge$edge_id, transform$edge_id), , drop = FALSE
    ]
    if (anyNA(transform$edge_id) ||
        any(!is.finite(transform$center)) ||
        any(!is.finite(transform$scale) | transform$scale <= 0)) {
      stop("ConditionGRNFit predictor transform is incomplete.",
           call. = FALSE)
    }
    beta <- as.matrix(fit$beta)
    contrast <- as.matrix(fit$contrast)
    mask <- as.matrix(fit$eligibility_mask)
    if (!identical(dim(beta), c(nrow(edge), length(fit$condition_levels))) ||
        !identical(dim(contrast), dim(beta)) ||
        !identical(dim(mask), dim(beta)) ||
        !identical(colnames(beta), fit$condition_levels) ||
        !identical(colnames(contrast), fit$condition_levels) ||
        !identical(colnames(mask), fit$condition_levels) ||
        !identical(rownames(beta), edge$edge_id) ||
        !identical(rownames(contrast), edge$edge_id) ||
        !identical(rownames(mask), edge$edge_id) ||
        any(!is.finite(beta)) || any(!is.finite(contrast)) ||
        !is.logical(mask) || anyNA(mask) ||
        !fit$reference_condition %in% colnames(beta)) {
      stop("ConditionGRNFit coefficient matrices are misaligned.",
           call. = FALSE)
    }
    beta_reference <- beta[, fit$reference_condition]
    expected_contrast <- .rc_reference_contrast(
      beta, fit$reference_condition
    )
    if (!isTRUE(all.equal(
      unname(contrast), unname(expected_contrast), tolerance = 1e-10
    ))) {
      stop("ConditionGRNFit reference contrasts are inconsistent.",
           call. = FALSE)
    }
    target_rsq <- as.matrix(fit$target_rsq)
    if (!all(c("target", "center", "scale") %in% colnames(response))) {
      stop("ConditionGRNFit target transforms are incomplete.",
           call. = FALSE)
    }
    response$target <- toupper(trimws(as.character(response$target)))
    if (anyDuplicated(response$target) ||
        any(!is.finite(response$center)) ||
        any(!is.finite(response$scale) | response$scale <= 0) ||
        is.null(rownames(target_rsq)) ||
        !identical(colnames(target_rsq), fit$condition_levels)) {
      stop("ConditionGRNFit target transforms or diagnostics are invalid.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      tab <- edge
      tab$tf <- toupper(trimws(as.character(tab$tf)))
      tab$target <- toupper(trimws(as.character(tab$target)))
      tab$region <- trimws(as.character(tab$region))
      tab$condition_estimate <- as.numeric(beta[, condition])
      tab$reference_estimate <- as.numeric(beta_reference)
      tab$condition_effect <- as.numeric(contrast[, condition])
      tab$estimate <- tab$condition_estimate
      tab$corr <- NA_real_
      tab$eligible_in_condition <- as.logical(mask[, condition])
      tab$reference_condition <- fit$reference_condition
      tab$predictor_center <- transform$center
      tab$predictor_scale <- transform$scale
      response_index <- match(tab$target, response$target)
      rsq_index <- match(tab$target, toupper(rownames(target_rsq)))
      if (anyNA(response_index) || anyNA(rsq_index)) {
        stop(
          "ConditionGRNFit edge targets do not align with target transforms ",
          "and diagnostics.", call. = FALSE
        )
      }
      tab$response_center <- response$center[response_index]
      tab$response_scale <- response$scale[response_index]
      tab$rsq <- target_rsq[rsq_index, condition]
      tab[[condition_col]] <- condition
      tab[[celltype_col]] <- fit$cell_type
      group_values <- tab[1L, c(condition_col, celltype_col), drop = FALSE]
      tab$group_id <- rc_make_stratum_id(
        group_values, c(condition_col, celltype_col)
      )
      tab$fit_engine <- fit$fit_engine
      tab$coefficient_scale <- fit$coefficient_scale
      tab <- tab[, c(
        "group_id", condition_col, celltype_col,
        setdiff(colnames(tab), c("group_id", condition_col, celltype_col))
      ), drop = FALSE]
      condition_rows[[row_index]] <- tab
      row_index <- row_index + 1L
    }
    summary <- edge
    summary$estimate <- rowMeans(beta)
    summary$corr <- NA_real_
    summary[[celltype_col]] <- fit$cell_type
    summary$reference_condition <- fit$reference_condition
    summary$summary_only <- TRUE
    universal_rows[[length(universal_rows) + 1L]] <- summary
  }
  all_edges <- do.call(rbind, condition_rows)
  rownames(all_edges) <- NULL
  condition_active <- all_edges[
    is.finite(all_edges$condition_estimate) &
      abs(all_edges$condition_estimate) >= active_tol &
      is.finite(all_edges$rsq) & all_edges$rsq >= min_model_rsq,
    , drop = FALSE
  ]
  effect_all <- all_edges
  effect_all$estimate <- effect_all$condition_effect
  effect_active <- effect_all[
    is.finite(effect_all$condition_effect) &
      abs(effect_all$condition_effect) >= active_tol &
      is.finite(effect_all$rsq) & effect_all$rsq >= min_model_rsq,
    , drop = FALSE
  ]
  params <- methods::slot(methods::slot(grn_object, "grn"), "params")
  list(
    network_index = params$condition_network_index %||% data.frame(),
    fit_diagnostics = params$condition_fit_diagnostics %||% data.frame(),
    fit_contracts = fits,
    universal = do.call(rbind, universal_rows),
    condition_all = all_edges,
    condition_active = condition_active,
    condition_effect_all = effect_all,
    condition_effect_active = effect_active,
    active_tol = active_tol
  )
}

.rc_run_condition_single_cell_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      method = "shared_design_independent",
      candidate_screen = "condition_union",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 1,
      condition_weight = "equal",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
      parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = FALSE,
    save_pando_objects = TRUE,
    BPPARAM = NULL,
    on_group_error = c("record", "stop"),
    species = c("auto", "human", "mouse")) {
  on_group_error <- match.arg(on_group_error)
  species <- .rc_infer_gem_species(gem, species)
  if (!is.numeric(min_cells) || length(min_cells) != 1L ||
      !is.finite(min_cells) || min_cells < 1 ||
      abs(min_cells - round(min_cells)) > sqrt(.Machine$double.eps)) {
    stop("`min_cells` must be one positive integer.", call. = FALSE)
  }
  min_cells <- as.integer(min_cells)
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj
  )
  if (!is.list(pando_initiate_args) || !is.list(pando_motif_args) ||
      !is.list(pando_infer_args)) {
    stop("Pando initiate, motif, and inference arguments must be lists.",
         call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop(
      "Install 1667857557/Pando_regcompass before condition-aware GRN inference.",
      call. = FALSE
    )
  }
  if (!exists("infer_condition_grn", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "Installed Pando lacks `infer_condition_grn`; install the v1.2.0 companion PR.",
      call. = FALSE
    )
  }
  if (!exists("condition_grn_fit", envir = asNamespace("Pando"),
              inherits = FALSE)) {
    stop(
      "Installed Pando lacks the ConditionGRNFit v2 accessor.",
      call. = FALSE
    )
  }
  pando_infer_args <- utils::modifyList(
    list(
      method = "shared_design_independent",
      candidate_screen = "condition_union",
      tf_cor = 0.1,
      peak_cor = 0.01,
      alpha = 0.5,
      condition_mix = 1,
      condition_weight = "equal",
      nlambda = 50L,
      nfolds = 5L,
      lambda_selection = "lambda.1se",
      scale = TRUE,
      parallel = FALSE
    ),
    pando_infer_args
  )
  required_infer <- list(
    method = "shared_design_independent",
    condition_mix = 1,
    condition_weight = "equal",
    scale = TRUE
  )
  invalid_infer <- names(required_infer)[!vapply(
    names(required_infer),
    function(name) {
      if (identical(name, "condition_mix")) {
        return(isTRUE(all.equal(
          suppressWarnings(as.numeric(pando_infer_args[[name]])), 1
        )))
      }
      if (identical(name, "scale")) {
        return(isTRUE(pando_infer_args[[name]]))
      }
      identical(
        as.character(pando_infer_args[[name]]),
        as.character(required_infer[[name]])
      )
    },
    logical(1)
  )]
  if (length(invalid_infer)) {
    stop(
      "RegCompass condition comparability requires Pando settings: ",
      paste(
        paste0(names(required_infer), "=", unlist(required_infer)),
        collapse = ", "
      ),
      ". Incompatible fields: ", paste(invalid_infer, collapse = ", "),
      call. = FALSE
    )
  }
  if (isTRUE(require_padj)) {
    warning(
      "`require_padj = TRUE` is ignored by the regularized condition solver; ",
      "active coefficients and target-level R-squared are used.",
      call. = FALSE
    )
  }
  pando_install <- .rc_validate_pando_repository()
  motif_policy <- "user_supplied"
  if (is.null(pfm)) {
    pfm <- .rc_default_pando_motifs()
    motif_policy <- "Pando::motifs"
  }
  if (!length(pfm)) stop("`pfm` must be non-empty.", call. = FALSE)
  group_cols <- c(condition_col, celltype_col)
  missing <- setdiff(group_cols, colnames(object@meta.data))
  if (length(missing)) {
    stop("Missing metadata columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop(
      "Pando condition inference requires cell-type-shared ATAC TF-IDF across conditions.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "pando_objects"),
             recursive = TRUE, showWarnings = FALSE)

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

  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Pando unified condition GRN"
  )
  object <- filtered$object
  init_defaults <- list(
    object = object,
    peak_assay = atac_assay,
    rna_assay = rna_assay
  )
  init_defaults[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init_defaults, pando_initiate_args))
  motif_defaults <- list(object = grn, pfm = pfm, genome = genome)
  motif_defaults[names(pando_motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif_defaults, pando_motif_args))

  infer_defaults <- list(
    object = grn,
    cell_type_col = celltype_col,
    condition_col = condition_col,
    genes = target_genes,
    network_name = "regcompass_condition_grn",
    min_cells_per_condition = min_cells,
    on_small_condition = if (identical(on_group_error, "stop")) "error" else
      "skip_cell_type",
    BPPARAM = if (identical(BPPARAM, FALSE)) NULL else BPPARAM
  )
  infer_defaults[names(pando_infer_args)] <- NULL
  grn <- do.call(
    Pando::infer_condition_grn,
    c(infer_defaults, pando_infer_args)
  )

  extracted <- .rc_extract_condition_grn_contract(
    grn,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq
  )
  meta <- object@meta.data
  expected <- unique(meta[, group_cols, drop = FALSE])
  expected$n_cells <- vapply(seq_len(nrow(expected)), function(i) {
    sum(
      as.character(meta[[condition_col]]) ==
        as.character(expected[[condition_col]][[i]]) &
      as.character(meta[[celltype_col]]) ==
        as.character(expected[[celltype_col]][[i]])
    )
  }, integer(1))
  expected$group_id <- rc_make_stratum_id(expected, group_cols)
  condition_index <- extracted$network_index[
    extracted$network_index$network_level == "condition", , drop = FALSE
  ]
  index_key <- paste(
    as.character(condition_index$condition),
    as.character(condition_index$cell_type),
    sep = "\001"
  )
  expected_key <- paste(
    as.character(expected[[condition_col]]),
    as.character(expected[[celltype_col]]),
    sep = "\001"
  )
  index_match <- match(expected_key, index_key)
  status <- expected
  status$n_target_genes <- length(target_genes)
  status$n_atac_peaks_input <- filtered$diagnostics$n_input_peaks
  status$n_zero_count_peaks_excluded <-
    filtered$diagnostics$n_zero_count_peaks_excluded
  status$n_atac_peaks_used <- filtered$diagnostics$n_retained_peaks
  status$status <- ifelse(is.na(index_match), "failed_missing_condition_network", "ok")
  status$n_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_all$group_id == id)
  }, integer(1))
  status$n_significant_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_active$group_id == id)
  }, integer(1))
  status$n_condition_effect_edges <- vapply(status$group_id, function(id) {
    sum(extracted$condition_effect_active$group_id == id)
  }, integer(1))
  status$error_class <- ifelse(status$status == "ok", NA_character_,
                               "missing_condition_network")
  status$error_message <- ifelse(
    status$status == "ok", NA_character_,
    "Pando did not emit a condition network for this condition-by-cell-type group."
  )
  status <- status[, c(
    "group_id", group_cols,
    setdiff(colnames(status), c("group_id", group_cols))
  ), drop = FALSE]
  if (any(status$status != "ok")) {
    stop(
      "Unified Pando condition GRN coverage is incomplete: ",
      paste(status$group_id[status$status != "ok"], collapse = "; "),
      call. = FALSE
    )
  }

  .rc_mm_write_tsv_gz(
    status, file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_active,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_effect_all,
    file.path(outdir, "pando_tf_peak_gene_condition_effect_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$condition_effect_active,
    file.path(outdir, "pando_tf_peak_gene_condition_effect_active.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$universal,
    file.path(outdir, "pando_tf_peak_gene_universal.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    extracted$network_index,
    file.path(outdir, "pando_condition_network_index.tsv.gz")
  )
  if (is.data.frame(extracted$fit_diagnostics) &&
      nrow(extracted$fit_diagnostics)) {
    .rc_mm_write_tsv_gz(
      extracted$fit_diagnostics,
      file.path(outdir, "pando_condition_fit_diagnostics.tsv.gz")
    )
  }
  predictor_transforms <- do.call(rbind, lapply(
    extracted$fit_contracts,
    function(fit) {
      out <- merge(
        fit$edge_table,
        fit$predictor_transform,
        by = "edge_id",
        all.x = TRUE,
        sort = FALSE
      )
      out$cell_type <- fit$cell_type
      out$reference_condition <- fit$reference_condition
      out$coefficient_scale <- fit$coefficient_scale
      out
    }
  ))
  .rc_mm_write_tsv_gz(
    predictor_transforms,
    file.path(outdir, "pando_edge_predictor_transforms.tsv.gz")
  )
  saveRDS(
    extracted$fit_contracts,
    file.path(outdir, "pando_condition_grn_fits.rds")
  )
  if (isTRUE(save_pando_objects)) {
    saveRDS(grn, file.path(outdir, "pando_objects",
                           "condition_grn_fit_v2.rds"))
  }

  tf_metabolic_target_overlap <- intersect(
    unique(toupper(extracted$condition_effect_active$tf)),
    unique(toupper(target_genes))
  )
  answer <- list(
    schema_version = "regcompass_condition_grn_fit_v2",
    pando_installed_version = pando_install$version,
    pando_installation = pando_install,
    target_metabolic_genes = target_genes,
    sample_status = status,
    pando_network_index = extracted$network_index,
    pando_fit_diagnostics = extracted$fit_diagnostics,
    condition_grn_fits = extracted$fit_contracts,
    tf_peak_gene_universal = extracted$universal,
    tf_peak_gene_condition_all = extracted$condition_all,
    tf_peak_gene_condition = extracted$condition_active,
    tf_peak_gene_condition_effect_all = extracted$condition_effect_all,
    tf_peak_gene_condition_effect = extracted$condition_effect_active,
    tf_metabolic_target_overlap = tf_metabolic_target_overlap,
    tf_peak_gene_all = extracted$condition_all,
    tf_peak_gene_significant = extracted$condition_active,
    normalization_policy = list(
      rna = "global single-cell normalized RNA shared by all conditions",
      atac = "cell-type-shared TF-IDF across conditions",
      grn_fit = paste(
        "one Pando fit contract per cell type with one TF-peak-target",
        "edge dictionary and independently estimated condition coefficients"
      ),
      universal_coefficient = "visualization-only row mean; never used as a contrast baseline",
      reference_condition = unique(vapply(
        extracted$fit_contracts, `[[`, character(1), "reference_condition"
      )),
      condition_effect = "condition coefficient minus explicit reference-condition coefficient",
      coefficient_scale =
        "pooled cell-type TF-by-ATAC edge and target standardization",
      core_reaction_evidence =
        "active condition-level TF-peak-target coefficients",
      penalty_regulatory_evidence = paste(
        "condition-minus-reference coefficients projected through metacell",
        "TF RNA by peak ATAC activity using Pando's stored edge transform"
      ),
      pando_motifs = motif_policy,
      pando_regions = region_policy,
      pando_padj = paste(
        "not applicable to the regularized condition solver;",
        "active coefficients and target-level R-squared are used"
      ),
      pando_peak_cor = pando_infer_args$peak_cor,
      legacy_padj_threshold = padj_threshold,
      legacy_require_padj = require_padj,
      active_tolerance = extracted$active_tol,
      min_model_rsq = min_model_rsq
    ),
    group_cols = group_cols
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_tf_peak_interaction <- function(tf_activity, peak_activity) {
  tf_activity <- as.matrix(tf_activity)
  peak_activity <- as.matrix(peak_activity)
  if (!identical(dim(tf_activity), dim(peak_activity))) {
    stop("TF RNA and peak ATAC edge matrices must have identical dimensions.",
         call. = FALSE)
  }
  if (!identical(colnames(tf_activity), colnames(peak_activity))) {
    stop("TF RNA and peak ATAC matrices must contain identically ordered units.",
         call. = FALSE)
  }
  tf_activity * peak_activity
}

.rc_project_condition_edges <- function(
    standardized_edge_activity, coefficient_contrast, target_rsq) {
  standardized_edge_activity <- as.matrix(standardized_edge_activity)
  coefficient_contrast <- suppressWarnings(as.numeric(coefficient_contrast))
  if (nrow(standardized_edge_activity) != length(coefficient_contrast) ||
      any(!is.finite(standardized_edge_activity)) ||
      any(!is.finite(coefficient_contrast))) {
    stop(
      "Standardized edge activity and coefficient contrast must align and be finite.",
      call. = FALSE
    )
  }
  target_rsq <- suppressWarnings(as.numeric(target_rsq))
  target_rsq <- target_rsq[is.finite(target_rsq)]
  reliability <- if (length(target_rsq)) {
    sqrt(min(max(stats::median(target_rsq), 0), 1))
  } else {
    0
  }
  reliability * as.numeric(
    crossprod(coefficient_contrast, standardized_edge_activity)
  )
}

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    target_genes = NULL) {
  if (!is.data.frame(significant_edges)) {
    stop("`significant_edges` must be a data.frame.", call. = FALSE)
  }
  required_edges <- c(
    "edge_id", "target", "region", "tf", "estimate",
    "predictor_center", "predictor_scale", "response_scale",
    condition_col, celltype_col
  )
  missing_edges <- setdiff(required_edges, colnames(significant_edges))
  if (length(missing_edges)) {
    stop(
      "Condition-effect edge table is missing columns: ",
      paste(missing_edges, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.data.frame(unit_meta) ||
      !all(c("pool_id", condition_col, celltype_col) %in% colnames(unit_meta))) {
    stop("`unit_meta` is incomplete for condition-pooled regulatory scoring.",
         call. = FALSE)
  }
  units <- colnames(object)
  unit_meta <- unit_meta[
    match(units, as.character(unit_meta$pool_id)), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata do not align to the scoring object.", call. = FALSE)
  }
  genes <- unique(tolower(trimws(as.character(target_genes))))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) {
    genes <- unique(tolower(trimws(as.character(significant_edges$target))))
  }
  modifier <- matrix(
    0,
    nrow = length(genes),
    ncol = length(units),
    dimnames = list(genes, units)
  )
  raw_projection <- modifier
  attr(modifier, "reliability_policy") <- paste(
    "finite condition-network R-squared is mapped through sqrt(clamp(median R2,0,1));",
    "targets without finite R-squared receive regulatory reliability zero"
  )
  attr(modifier, "regulatory_predictor") <- paste(
    "tanh of the unnormalised model-space projection:",
    "sum_edge standardized(TF_RNA x peak_ATAC) * (beta_condition-beta_reference)"
  )
  attr(modifier, "raw_model_projection") <- raw_projection
  attr(modifier, "projection_formula") <- paste(
    "sqrt(clamp(target_R2,0,1)) * sum_edge",
    "[(TF_RNA*peak_ATAC-center_edge)/scale_edge] *",
    "(beta_condition-beta_reference); no L1 edge-weight normalisation"
  )
  if (!nrow(significant_edges) || !length(genes)) return(modifier)

  edges <- significant_edges
  rsq <- if ("rsq" %in% colnames(edges)) {
    suppressWarnings(as.numeric(edges$rsq))
  } else {
    rep(NA_real_, nrow(edges))
  }
  edges <- edges[is.finite(rsq), , drop = FALSE]
  if (!nrow(edges)) return(modifier)
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  edges$predictor_center <- suppressWarnings(
    as.numeric(edges$predictor_center)
  )
  edges$predictor_scale <- suppressWarnings(
    as.numeric(edges$predictor_scale)
  )
  edges$response_scale <- suppressWarnings(
    as.numeric(edges$response_scale)
  )
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & edges$estimate != 0 &
      is.finite(edges$predictor_center) &
      is.finite(edges$predictor_scale) & edges$predictor_scale > 0 &
      is.finite(edges$response_scale) & edges$response_scale > 0,
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  rna <- .rc_pando_assay_data(object, rna_assay)
  atac <- .rc_pando_assay_data(object, atac_assay)
  if (!setequal(colnames(rna), units) || !setequal(colnames(atac), units)) {
    stop("Normalized RNA/ATAC data and metacell units do not match.",
         call. = FALSE)
  }
  rna <- rna[, units, drop = FALSE]
  atac <- atac[, units, drop = FALSE]
  tf_keys <- toupper(trimws(rownames(rna)))
  tf_keep <- !is.na(tf_keys) & nzchar(tf_keys) & !duplicated(tf_keys)
  tf_lookup <- stats::setNames(rownames(rna)[tf_keep], tf_keys[tf_keep])
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  peak_keep <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[peak_keep], peak_keys[peak_keep])
  edges$.tf_id <- unname(tf_lookup[edges$tf])
  edges$.peak_id <- unname(
    peak_lookup[toupper(.rc_pando_region_key(edges$region))]
  )
  edges <- edges[
    !is.na(edges$.tf_id) & nzchar(edges$.tf_id) &
      !is.na(edges$.peak_id) & nzchar(edges$.peak_id),
    , drop = FALSE
  ]
  if (!nrow(edges)) return(modifier)

  group_key_edges <- paste(
    as.character(edges[[condition_col]]),
    as.character(edges[[celltype_col]]),
    sep = "\001"
  )
  group_key_units <- paste(
    as.character(unit_meta[[condition_col]]),
    as.character(unit_meta[[celltype_col]]),
    sep = "\001"
  )
  for (group_key in unique(group_key_edges)) {
    group_edges <- edges[group_key_edges == group_key, , drop = FALSE]
    group_units <- units[group_key_units == group_key]
    if (!nrow(group_edges) || !length(group_units)) next
    for (target in unique(group_edges$target)) {
      selected <- group_edges[group_edges$target == target, , drop = FALSE]
      gene_id <- tolower(target)
      if (!gene_id %in% rownames(modifier) || !nrow(selected)) next
      tf_activity <- as.matrix(
        rna[selected$.tf_id, units, drop = FALSE]
      )
      peak_activity <- as.matrix(
        atac[selected$.peak_id, units, drop = FALSE]
      )
      edge_activity <- .rc_tf_peak_interaction(tf_activity, peak_activity)
      edge_model <- sweep(
        edge_activity, 1L, selected$predictor_center, "-"
      )
      edge_model <- sweep(
        edge_model, 1L, selected$predictor_scale, "/"
      )
      edge_model <- edge_model[, group_units, drop = FALSE]
      value_raw <- .rc_project_condition_edges(
        standardized_edge_activity = edge_model,
        coefficient_contrast = selected$estimate,
        target_rsq = selected$rsq
      )
      raw_projection[gene_id, group_units] <- value_raw
      modifier[gene_id, group_units] <- tanh(value_raw)
    }
  }
  attr(modifier, "tf_target_overlap") <- intersect(
    unique(edges$tf), unique(edges$target)
  )
  attr(modifier, "raw_model_projection") <- raw_projection
  modifier
}

.rc_integrate_regulatory_support <- function(
    rna_support, regulatory_modifier, alpha = 1) {
  rna_support <- as.matrix(rna_support)
  regulatory_modifier <- as.matrix(regulatory_modifier)
  if (!identical(dim(rna_support), dim(regulatory_modifier)) ||
      !identical(dimnames(rna_support), dimnames(regulatory_modifier))) {
    stop("RNA support and regulatory modifier matrices must align exactly.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0) {
    stop("`alpha` must be one finite non-negative number.", call. = FALSE)
  }
  C <- pmin(pmax(rna_support, 0), 1)
  R <- pmin(pmax(regulatory_modifier, -1), 1)
  multiplier <- 2^(alpha * R)
  numerator <- C * multiplier
  denominator <- 1 - C + numerator
  out <- numerator / denominator
  out[C <= 0] <- 0
  out[C >= 1] <- 1
  out[!is.finite(out)] <- NA_real_
  dimnames(out) <- dimnames(C)
  attr(out, "integration_formula") <- paste(
    "C_multiome = C_RNA * 2^(alpha * R_condition_TFxATAC) /",
    "(1 - C_RNA + C_RNA * 2^(alpha * R_condition_TFxATAC))"
  )
  attr(out, "score_semantics") <- paste(
    "zero-preserving bounded target-gene RNA support with a signed",
    "condition-effect TF-by-ATAC modifier on the support log-odds scale"
  )
  out
}

.rc_build_condition_pooled_layer1 <- function(
    metacell_object, meta_modules, gem, metacell_meta,
    sample_col = "sample_id", condition_col = "condition",
    celltype_col = "cell_type", rna_assay = "RNA", atac_assay = "ATAC",
    regulatory_alpha = 1,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = NULL) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  gpr_and_method <- match.arg(gpr_and_method)
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  full_library_size <- Matrix::colSums(counts)
  keep <- tolower(rownames(counts)) %in% gpr_genes
  rna_counts <- counts[keep, , drop = FALSE]
  rna_logcpm <- .rc_metacell_logcpm(
    rna_counts,
    library_size = full_library_size[colnames(rna_counts)]
  )
  rownames(rna_logcpm) <- tolower(rownames(rna_logcpm))
  if (anyDuplicated(rownames(rna_logcpm))) {
    stop("Duplicated GPR gene identifiers after case normalization.",
         call. = FALSE)
  }

  unit_meta <- metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("Pooled metacell metadata lack metacell_id/pool_id.", call. = FALSE)
  }
  unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  unit_meta$unit_id <- unit_meta$pool_id
  unit_meta[[sample_col]] <- paste0(
    as.character(unit_meta[[condition_col]]), "__pooled"
  )
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), unit_meta$pool_id), , drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Pooled metacell metadata do not align with RNA counts.",
         call. = FALSE)
  }

  gene_rna_support <- rc_gene_score(
    rna_logcpm,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  regulatory_edges <- meta_modules$tf_peak_gene_condition_effect
  if (!is.data.frame(regulatory_edges)) {
    stop(
      "Stage 3 lacks condition-effect TF-peak-target edges. Rerun GRN and ",
      "meta-module stages with the unified condition-aware Pando workflow.",
      call. = FALSE
    )
  }
  modifier <- .rc_condition_gene_regulatory_modifier(
    significant_edges = regulatory_edges,
    object = metacell_object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    target_genes = rownames(gene_rna_support)
  )
  raw_projection <- attr(modifier, "raw_model_projection")
  if (is.null(raw_projection)) {
    raw_projection <- matrix(
      0, nrow = nrow(modifier), ncol = ncol(modifier),
      dimnames = dimnames(modifier)
    )
  }
  modifier <- modifier[
    rownames(gene_rna_support),
    colnames(gene_rna_support),
    drop = FALSE
  ]
  raw_projection <- raw_projection[
    rownames(gene_rna_support),
    colnames(gene_rna_support),
    drop = FALSE
  ]
  gene_multiome_support <- .rc_integrate_regulatory_support(
    gene_rna_support,
    modifier,
    alpha = regulatory_alpha
  )
  reaction_expression <- rc_reaction_capacity(
    parsed,
    gene_multiome_support,
    promiscuity_mode = "none",
    and_method = gpr_and_method,
    or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )

  list(
    schema_version = "regcompass_condition_grn_layer1_v2",
    reaction_expression = reaction_expression,
    rna_metacell_logcpm = rna_logcpm,
    gene_support_rna = gene_rna_support,
    gene_regulatory_modifier = modifier,
    gene_regulatory_model_projection = raw_projection,
    gene_support_multiome = gene_multiome_support,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, rownames(rna_logcpm)),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "condition_only_metacell_with_posthoc_celltype",
    capacity_params = list(
      regulatory_alpha = regulatory_alpha,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode =
        "fixed_pando_transform_x_reference_coefficient_contrast",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel,
      bpparam_class = if (is.null(BPPARAM)) {
        "auto_or_sequential"
      } else if (identical(BPPARAM, FALSE)) {
        "sequential"
      } else {
        class(BPPARAM)[[1L]]
      }
    ),
    evidence_formula = paste(
      "condition GRN active edges define metabolic targets/core reactions;",
      "Pando reference contrasts project the identically transformed TF RNA x peak ATAC predictor;",
      "edge effects are summed without L1 normalisation and bounded once with tanh;",
      "the signed modifier updates target-gene RNA support on the log-odds scale;",
      paste0("COMPASS-compatible ", gpr_and_method, " GPR-AND"),
      "and additive isozyme OR"
    ),
    evidence_inputs = c(
      "target_metabolic_gene_RNA",
      "regulatory_TF_RNA",
      "regulatory_peak_ATAC",
      "condition_level_Pando_coefficient",
      "reference_condition_Pando_coefficient",
      "condition_minus_reference_coefficient",
      "Pando_pooled_edge_center_and_scale"
    )
  )
}
