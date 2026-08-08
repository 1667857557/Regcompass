# RegCompass-owned condition x cell-type scheduling over public Pando primitives.
#
# The mathematical contract remains Pando's canonical two-stage procedure:
# global + per-condition candidate discovery -> exact edge union within each
# cell type -> fixed-dictionary Gaussian identity GLM in every condition.
# RegCompass changes only task scheduling, worker lifetime and result assembly.

.rc_condition_network_label <- function(x) {
  value <- gsub("[^[:alnum:]_.-]+", "_", as.character(x))
  value[!nzchar(value)] <- "unnamed"
  value
}

.rc_condition_parallel_plan <- function(
    metadata, condition_types, condition_col, celltype_col, min_cells) {
  condition_types <- unique(as.character(condition_types))
  if (!length(condition_types)) {
    stop("No condition-GRN cell type was supplied.", call. = FALSE)
  }
  plans <- vector("list", length(condition_types))
  names(plans) <- condition_types
  for (type in condition_types) {
    type_cells <- rownames(metadata)[
      as.character(metadata[[celltype_col]]) == type
    ]
    levels <- unique(as.character(metadata[type_cells, condition_col]))
    counts <- vapply(levels, function(level) {
      sum(as.character(metadata[type_cells, condition_col]) == level)
    }, integer(1))
    undersized <- levels[counts < as.integer(min_cells)]
    if (length(undersized)) {
      detail <- paste0(undersized, "=", counts[undersized], collapse = ", ")
      stop(
        "Cell type `", type, "` has condition(s) below min_cells: ",
        detail, call. = FALSE
      )
    }
    if (length(levels) < 2L) {
      stop(
        "Condition-GRN cell type `", type,
        "` must retain at least two conditions.", call. = FALSE
      )
    }
    cells_by_condition <- stats::setNames(lapply(levels, function(level) {
      type_cells[as.character(metadata[type_cells, condition_col]) == level]
    }), levels)
    plans[[type]] <- list(
      cell_type = type,
      conditions = levels,
      cells_by_condition = cells_by_condition,
      global_cells = unlist(cells_by_condition, use.names = FALSE)
    )
  }
  plans
}

.rc_condition_discovery_task <- function(
    task, grn, target_genes, pando_infer_args) {
  thread_state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(thread_state), add = TRUE)
  edge <- Pando::discover_grn_edges(
    object = grn,
    genes = target_genes,
    cells = task$cells,
    source_label = task$source_label,
    source_type = task$source_type,
    tf_cor = pando_infer_args$tf_cor,
    peak_cor = pando_infer_args$peak_cor,
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    parallel = FALSE,
    verbose = FALSE
  )
  list(
    cell_type = task$cell_type,
    condition = task$condition,
    source_type = task$source_type,
    edge = edge
  )
}

.rc_condition_fit_task <- function(
    task, grn, pando_infer_args) {
  thread_state <- .rc_set_internal_single_thread()
  on.exit(.rc_restore_internal_threads(thread_state), add = TRUE)
  fitted <- Pando::fit_grn_from_edges(
    object = grn,
    edge_dictionary = task$dictionary,
    cells = task$cells,
    condition_label = task$condition,
    network_name = task$network_name,
    adjust_method = "BH",
    padj_threshold = 0.05,
    rank_action = pando_infer_args$rank_action,
    min_residual_df = pando_infer_args$min_residual_df,
    rna_layer = pando_infer_args$rna_layer %||% "data",
    peak_layer = pando_infer_args$peak_layer %||% "data",
    peak_value_type = pando_infer_args$peak_value_type %||% "normalized",
    parallel = FALSE,
    overwrite = TRUE,
    verbose = FALSE
  )
  network <- Pando::GetNetwork(fitted, network = task$network_name)
  list(
    cell_type = task$cell_type,
    condition = task$condition,
    network_name = task$network_name,
    network = network,
    coefficients = as.data.frame(stats::coef(network), stringsAsFactors = FALSE),
    fit = as.data.frame(Pando::gof(network), stringsAsFactors = FALSE)
  )
}

.rc_condition_attach_parallel_contract <- function(
    grn, plans, dictionaries, fit_results, condition_col, celltype_col,
    pando_infer_args, parallel_plan) {
  regulatory <- methods::slot(grn, "grn")
  networks <- methods::slot(regulatory, "networks")
  params <- methods::slot(regulatory, "params")
  index <- list()

  for (result in fit_results) {
    if (result$network_name %in% names(networks)) {
      stop(
        "Duplicated condition network during RegCompass merge: ",
        result$network_name, call. = FALSE
      )
    }
    networks[[result$network_name]] <- result$network
    dictionary <- dictionaries[[result$cell_type]]
    index[[length(index) + 1L]] <- data.frame(
      cell_type = result$cell_type,
      condition = result$condition,
      network_name = result$network_name,
      n_cells = length(
        plans[[result$cell_type]]$cells_by_condition[[result$condition]]
      ),
      n_dictionary_edges = nrow(dictionary),
      n_significant_edges = sum(
        result$coefficients$significant %in% TRUE, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
  methods::slot(regulatory, "networks") <- networks
  if (length(fit_results)) {
    methods::slot(regulatory, "active_network") <-
      fit_results[[length(fit_results)]]$network_name
  }

  fits <- list()
  for (type in names(plans)) {
    conditions <- plans[[type]]$conditions
    one <- fit_results[vapply(fit_results, function(x) {
      identical(x$cell_type, type)
    }, logical(1))]
    names(one) <- vapply(one, `[[`, character(1), "condition")
    one <- one[conditions]
    if (any(vapply(one, is.null, logical(1)))) {
      stop("Condition fixed-dictionary result set is incomplete.", call. = FALSE)
    }
    coefficient <- do.call(rbind, lapply(one, `[[`, "coefficients"))
    fit_table <- do.call(rbind, lapply(one, `[[`, "fit"))
    rownames(coefficient) <- rownames(fit_table) <- NULL
    dictionary <- dictionaries[[type]]
    network_names <- stats::setNames(
      vapply(one, `[[`, character(1), "network_name"), conditions
    )
    fit_contract <- list(
      schema_version = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
      fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
      coefficient_scale = "shared_preprocessed_input_units_unscaled",
      inference_scope = "conditional_on_selected_edge_dictionary",
      cell_type = type,
      condition_levels = conditions,
      condition_col = condition_col,
      cell_type_col = celltype_col,
      condition_cell_ids = plans[[type]]$cells_by_condition,
      edge_dictionary = dictionary,
      coefficients = coefficient,
      fit = fit_table,
      network_names = network_names,
      padj_threshold = 0.05,
      adjust_method = "BH",
      scale = FALSE,
      interaction = ":",
      projection_effect_column = "penalty_effect",
      projection_policy = "padj_significant_effects_only",
      target_genes = unique(as.character(dictionary$target)),
      rna_assay = Pando::Params(grn)$rna_assay,
      atac_assay = Pando::Params(grn)$peak_assay,
      rna_layer = attr(dictionary, "rna_layer", exact = TRUE) %||%
        (pando_infer_args$rna_layer %||% "data"),
      peak_layer = attr(dictionary, "peak_layer", exact = TRUE) %||%
        (pando_infer_args$peak_layer %||% "data"),
      peak_value_type = attr(
        dictionary, "peak_value_type", exact = TRUE
      ) %||% (pando_infer_args$peak_value_type %||% "normalized"),
      preprocessing_fingerprint = attr(
        dictionary, "preprocessing_fingerprint", exact = TRUE
      ),
      dictionary_preprocessing_provenance_verified = isTRUE(attr(
        dictionary, "preprocessing_provenance_verified", exact = TRUE
      ))
    )
    class(fit_contract) <- c("ConditionGRNFit", "list")
    invisible(.rc_require_pando_condition_grn_fit_schema(fit_contract))
    fits[[type]] <- fit_contract
  }

  params$analysis_mode <- "condition_grn"
  params$condition_col <- condition_col
  params$condition_levels <- unique(unlist(lapply(
    plans, `[[`, "conditions"), use.names = FALSE
  ))
  params$cell_type_col <- celltype_col
  params$condition_coefficients_calculated <- TRUE
  params$condition_grn_schema <- .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
  params$condition_grn_fits <- fits
  params$condition_network_index <- do.call(rbind, index)
  params$parallel_plan <- parallel_plan
  methods::slot(regulatory, "params") <- params
  methods::slot(grn, "grn") <- regulatory
  grn
}

.rc_fit_condition_grns_regcompass_parallel <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    cell_type = NULL, rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      tf_cor = 0.1, peak_cor = 0, adjust_method = "BH",
      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L
    ),
    save_pando_objects = TRUE, BPPARAM = NULL,
    progress_monitor = NULL,
    species = c("auto", "human", "mouse")) {
  species <- .rc_infer_gem_species(gem, species)
  rc_validate_gem(gem)
  .rc_validate_condition_celltype_metadata(
    object@meta.data, condition_col, celltype_col,
    require_multiple_conditions = TRUE
  )
  if (!is.list(pando_infer_args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  defaults <- list(
    tf_cor = 0.1, peak_cor = 0, adjust_method = "BH",
    padj_threshold = 0.05, rank_action = "mark", min_residual_df = 1L,
    rna_layer = "data", peak_layer = "data",
    peak_value_type = "normalized"
  )
  pando_infer_args <- utils::modifyList(defaults, pando_infer_args)
  if (!identical(toupper(as.character(pando_infer_args$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(pando_infer_args$padj_threshold), 0.05))) {
    stop("Canonical RegCompass condition effects require BH padj < 0.05.",
         call. = FALSE)
  }
  condition_types <- if (is.null(cell_type)) {
    unique(as.character(object@meta.data[[celltype_col]]))
  } else {
    unique(as.character(cell_type))
  }
  missing_types <- setdiff(
    condition_types, unique(as.character(object@meta.data[[celltype_col]]))
  )
  if (length(missing_types)) {
    stop(
      "Condition-GRN cell type(s) absent from Stage 1 object: ",
      paste(missing_types, collapse = ", "), call. = FALSE
    )
  }

  plans <- .rc_condition_parallel_plan(
    metadata = object@meta.data,
    condition_types = condition_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_cells = min_cells
  )
  if (is.null(pfm)) pfm <- .rc_default_pando_motifs()
  if (!"regions" %in% names(pando_initiate_args) ||
      is.null(pando_initiate_args$regions)) {
    pando_initiate_args$regions <- .rc_default_pando_regions(species)
  }
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  rna_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  target_upper <- intersect(toupper(rna_genes), toupper(metabolic_genes))
  target_genes <- rna_genes[toupper(rna_genes) %in% target_upper]
  if (!length(target_genes)) {
    stop("No overlap between RNA genes and GEM metabolic genes.", call. = FALSE)
  }
  filtered <- .rc_drop_zero_count_atac_features(
    object, atac_assay, "Pando condition GRNs"
  )
  object <- filtered$object
  init <- list(object = object, peak_assay = atac_assay, rna_assay = rna_assay)
  init[names(pando_initiate_args)] <- NULL
  grn <- do.call(Pando::initiate_grn, c(init, pando_initiate_args))
  pando_motif_args <- .rc_regcompass_motif_args(pando_motif_args)
  motif <- list(object = grn, pfm = pfm, genome = genome)
  motif[names(pando_motif_args)] <- NULL
  grn <- do.call(Pando::find_motifs, c(motif, pando_motif_args))

  discovery_tasks <- list()
  for (type in names(plans)) {
    plan <- plans[[type]]
    discovery_tasks[[length(discovery_tasks) + 1L]] <- list(
      cell_type = type,
      condition = NA_character_,
      source_label = "global",
      source_type = "global",
      cells = plan$global_cells
    )
    for (condition in plan$conditions) {
      discovery_tasks[[length(discovery_tasks) + 1L]] <- list(
        cell_type = type,
        condition = condition,
        source_label = condition,
        source_type = "condition",
        cells = plan$cells_by_condition[[condition]]
      )
    }
  }
  discovery_param <- if (length(discovery_tasks) > 1L) BPPARAM else FALSE
  .rc_step_monitor_event(
    progress_monitor, "condition_candidate_discovery",
    "dispatching global and condition x cell-type Pando candidate discovery",
    current = 6L,
    context = list(tasks = length(discovery_tasks), nested_parallel = FALSE)
  )
  discovery <- rc_parallel_lapply(
    discovery_tasks,
    .rc_condition_discovery_task,
    BPPARAM = discovery_param,
    grn = grn,
    target_genes = target_genes,
    pando_infer_args = pando_infer_args
  )
  .rc_step_monitor_event(
    progress_monitor, "condition_candidate_discovery_complete",
    "all candidate tasks completed; worker pool released",
    current = 7L,
    context = list(tasks = length(discovery_tasks))
  )

  dictionaries <- list()
  for (type in names(plans)) {
    global <- discovery[vapply(discovery, function(x) {
      identical(x$cell_type, type) && identical(x$source_type, "global")
    }, logical(1))]
    by_condition <- discovery[vapply(discovery, function(x) {
      identical(x$cell_type, type) && identical(x$source_type, "condition")
    }, logical(1))]
    if (length(global) != 1L ||
        length(by_condition) != length(plans[[type]]$conditions)) {
      stop("Condition candidate discovery returned an incomplete task set.",
           call. = FALSE)
    }
    condition_edges <- lapply(by_condition, `[[`, "edge")
    names(condition_edges) <- vapply(
      by_condition, `[[`, character(1), "condition"
    )
    condition_edges <- condition_edges[plans[[type]]$conditions]
    dictionaries[[type]] <- Pando::union_grn_edges(
      global_edges = global[[1L]]$edge,
      condition_edges = condition_edges
    )
  }
  .rc_step_monitor_event(
    progress_monitor, "condition_dictionary_barrier",
    "froze one exact edge dictionary per cell type",
    current = 8L,
    context = list(cell_types = length(dictionaries))
  )
  invisible(gc(verbose = FALSE, full = TRUE))

  fit_tasks <- list()
  for (type in names(plans)) {
    for (condition in plans[[type]]$conditions) {
      fit_tasks[[length(fit_tasks) + 1L]] <- list(
        cell_type = type,
        condition = condition,
        cells = plans[[type]]$cells_by_condition[[condition]],
        dictionary = dictionaries[[type]],
        network_name = paste0(
          "regcompass_condition_grn__",
          .rc_condition_network_label(type),
          "__condition__",
          .rc_condition_network_label(condition)
        )
      )
    }
  }
  fit_param <- if (length(fit_tasks) > 1L) BPPARAM else FALSE
  .rc_step_monitor_event(
    progress_monitor, "fixed_dictionary_fit",
    "dispatching condition x cell-type fixed-dictionary Pando GLMs",
    current = 9L,
    context = list(tasks = length(fit_tasks), nested_parallel = FALSE)
  )
  fit_results <- rc_parallel_lapply(
    fit_tasks,
    .rc_condition_fit_task,
    BPPARAM = fit_param,
    grn = grn,
    pando_infer_args = pando_infer_args
  )
  .rc_step_monitor_event(
    progress_monitor, "fixed_dictionary_fit_complete",
    "all condition x cell-type GLMs completed; worker pool released",
    current = 9L,
    context = list(tasks = length(fit_tasks))
  )

  parallel_plan <- list(
    scope = "condition_x_cell_type",
    candidate_discovery_tasks = length(discovery_tasks),
    fixed_dictionary_fit_tasks = length(fit_tasks),
    nested_target_parallel = FALSE,
    worker_budget = .rc_bpparam_worker_limit(BPPARAM, default = 1L),
    stage_barrier =
      "candidate_discovery_then_exact_union_then_fixed_dictionary_fit",
    scheduler_owner = "RegCompassR"
  )
  grn <- .rc_condition_attach_parallel_contract(
    grn = grn,
    plans = plans,
    dictionaries = dictionaries,
    fit_results = fit_results,
    condition_col = condition_col,
    celltype_col = celltype_col,
    pando_infer_args = pando_infer_args,
    parallel_plan = parallel_plan
  )

  extracted <- .rc_extract_condition_grn_contract(
    grn, condition_col, celltype_col
  )
  execution_summary <- .rc_pando_execution_summary(
    extracted$fit_diagnostics
  )
  execution_summary$parallel_plan <- parallel_plan
  meta <- object@meta.data
  status_rows <- list()
  for (fit in extracted$fit_contracts) {
    for (condition in fit$condition_levels) {
      cells <- fit$condition_cell_ids[[condition]]
      key_frame <- data.frame(
        condition_value = condition,
        celltype_value = fit$cell_type,
        stringsAsFactors = FALSE
      )
      names(key_frame) <- c(condition_col, celltype_col)
      id <- rc_make_stratum_id(
        key_frame, c(condition_col, celltype_col)
      )
      all_rows <- extracted$condition_all[
        extracted$condition_all[[condition_col]] == condition &
          extracted$condition_all[[celltype_col]] == fit$cell_type,
        , drop = FALSE
      ]
      active_rows <- extracted$condition_active[
        extracted$condition_active[[condition_col]] == condition &
          extracted$condition_active[[celltype_col]] == fit$cell_type,
        , drop = FALSE
      ]
      one <- data.frame(
        group_id = id,
        condition_value = condition,
        celltype_value = fit$cell_type,
        n_cells = length(cells),
        status = "ok",
        n_target_genes = length(unique(all_rows$target)),
        n_edges = nrow(all_rows),
        n_active_edges = nrow(active_rows),
        grn_evidence_role =
          "within_cell_type_common_dictionary_condition_glm",
        stringsAsFactors = FALSE
      )
      names(one)[names(one) == "condition_value"] <- condition_col
      names(one)[names(one) == "celltype_value"] <- celltype_col
      status_rows[[length(status_rows) + 1L]] <- one
    }
  }
  status <- do.call(rbind, status_rows)
  rownames(status) <- NULL
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(status, file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(
    extracted$condition_all,
    file.path(outdir, "pando_tf_peak_gene_condition_all.tsv.gz")
  )
  .rc_write_tsv_gz(
    extracted$condition_active,
    file.path(outdir, "pando_tf_peak_gene_condition_active.tsv.gz")
  )
  .rc_write_tsv_gz(
    extracted$universal,
    file.path(outdir, "pando_tf_peak_gene_universal.tsv.gz")
  )
  saveRDS(
    extracted$fit_contracts,
    file.path(outdir, "pando_condition_grn_fits.rds")
  )
  if (isTRUE(save_pando_objects)) {
    dir.create(
      file.path(outdir, "pando_objects"),
      recursive = TRUE, showWarnings = FALSE
    )
    saveRDS(
      grn,
      file.path(outdir, "pando_objects", "condition_grn_fit.rds")
    )
  }
  selected_cells <- unique(unlist(lapply(
    extracted$fit_contracts, function(fit) {
      unlist(fit$condition_cell_ids, use.names = FALSE)
    }
  ), use.names = FALSE))
  answer <- list(
    schema_version = "regcompass_condition_grn_common_dictionary_v1",
    analysis_mode = "condition_grn",
    condition_coefficients_calculated = TRUE,
    pando_fit_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    pando_installed_version = as.character(utils::packageVersion("Pando")),
    pando_grn_data = grn,
    paired_cell_ids = selected_cells,
    paired_cell_metadata = data.frame(
      cell_id = selected_cells,
      condition = as.character(meta[selected_cells, condition_col]),
      cell_type = as.character(meta[selected_cells, celltype_col]),
      stringsAsFactors = FALSE
    ),
    target_metabolic_genes = target_genes,
    condition_fit_status = status,
    pando_network_index = extracted$network_index,
    pando_fit_diagnostics = extracted$fit_diagnostics,
    pando_execution_summary = execution_summary,
    condition_grn_fits = extracted$fit_contracts,
    tf_peak_gene_universal = extracted$universal,
    tf_peak_gene_condition_all = extracted$condition_all,
    tf_peak_gene_condition = extracted$condition_active,
    tf_peak_gene_condition_effect_all = extracted$condition_all,
    tf_peak_gene_condition_effect = extracted$condition_active,
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF across conditions",
      grn_fit = paste(
        "global-plus-condition candidate discovery, exact edge union,",
        "fixed-dictionary Gaussian GLM"
      ),
      condition_effect = "unscaled fixed-dictionary condition coefficient",
      coefficient_contract =
        "same_exact_edge_dictionary_unscaled_gaussian_glm",
      significance = "estimable and BH adjusted P below 0.05",
      parallel_contract = parallel_plan,
      penalty_regulatory_evidence = paste(
        "paired-cell TF-by-ATAC projection using penalty_effect without",
        "effect-size or model-R2 gates"
      )
    ),
    group_cols = c(condition_col, celltype_col)
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
