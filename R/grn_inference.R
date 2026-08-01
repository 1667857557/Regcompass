.rc_require_normalized_assay <- function(object, assay, context) {
  value <- .rc_pando_assay_data(object, assay)
  expected <- as.character(colnames(object))
  observed <- as.character(colnames(value))
  if (anyNA(expected) || any(!nzchar(expected)) || anyDuplicated(expected) ||
      anyNA(observed) || any(!nzchar(observed)) || anyDuplicated(observed)) {
    stop(
      context,
      " normalized assay cell identifiers must be unique and non-empty.",
      call. = FALSE
    )
  }
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  if (length(missing) || length(extra)) {
    stop(
      context,
      " normalized assay contains different cells from the analysis object. ",
      "Missing: ", paste(utils::head(missing, 10L), collapse = ", "),
      "; extra: ", paste(utils::head(extra, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_run_condition_single_cell_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      method = "glm", tf_cor = 0.1, peak_cor = 0.01,
      adjust_method = "fdr", parallel = FALSE
    ),
    padj_threshold = 0.05,
    min_abs_estimate = 0,
    min_model_rsq = 0.1,
    require_padj = TRUE,
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
  if (!is.logical(save_pando_objects) || length(save_pando_objects) != 1L ||
      is.na(save_pando_objects)) {
    stop("`save_pando_objects` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(pando_initiate_args)) {
    stop("`pando_initiate_args` must be a list.", call. = FALSE)
  }
  if (!is.list(pando_motif_args)) {
    stop("`pando_motif_args` must be a list.", call. = FALSE)
  }
  .rc_assert_pando_nested_boundary(
    pando_motif_args,
    c("object", "pfm", "genome", "store_motif_positions"),
    "pando_motif_args"
  )
  if (!is.list(pando_infer_args)) {
    stop("`pando_infer_args` must be a list.", call. = FALSE)
  }
  pando_infer_args <- modifyList(
    list(
      method = "glm", tf_cor = 0.1, peak_cor = 0.01,
      adjust_method = "fdr", parallel = FALSE
    ),
    pando_infer_args
  )
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop(
      paste(
        "Install Pando from 1667857557/Pando_regcompass, either from GitHub",
        "or a local source tarball, before running GRN inference."
      ),
      call. = FALSE
    )
  }
  pando_install <- .rc_validate_pando_repository()
  motif_policy <- "user_supplied"
  if (is.null(pfm)) {
    pfm <- .rc_default_pando_motifs()
    motif_policy <- "Pando::motifs"
  }
  if (!length(pfm)) {
    stop("`pfm` must be a non-empty motif collection.", call. = FALSE)
  }
  group_cols <- c(condition_col, celltype_col)
  missing <- setdiff(group_cols, colnames(object@meta.data))
  if (length(missing)) {
    stop(
      "Missing metadata columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  .rc_require_normalized_assay(object, rna_assay, "RNA")
  .rc_require_normalized_assay(object, atac_assay, "ATAC")
  normalization <- object@misc$regcompass_atac_normalization %||% list()
  if (!identical(normalization$scope, "cell_type_across_conditions")) {
    stop(
      "Pando requires cell-type-shared ATAC TF-IDF across conditions.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(
    file.path(outdir, "pando_objects"),
    recursive = TRUE,
    showWarnings = FALSE
  )

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
  meta$.rc_pando_group_id <- rc_make_stratum_id(meta, group_cols)
  group_ids <- unique(as.character(meta$.rc_pando_group_id))

  run_one_group <- function(group_id) {
    cells <- rownames(meta)[as.character(meta$.rc_pando_group_id) == group_id]
    vals <- meta[
      match(cells[[1L]], rownames(meta)), group_cols, drop = FALSE
    ]
    status <- data.frame(
      group_id = group_id,
      condition = as.character(vals[[condition_col]][[1L]]),
      cell_type = as.character(vals[[celltype_col]][[1L]]),
      n_cells = length(cells),
      n_target_genes = length(target_genes),
      n_atac_peaks_input = NA_integer_,
      n_zero_count_peaks_excluded = NA_integer_,
      n_atac_peaks_used = NA_integer_,
      status = "pending",
      n_edges = 0L,
      n_significant_edges = 0L,
      error_class = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    names(status)[2:3] <- group_cols
    if (length(cells) < min_cells) {
      status$status <- "skipped_too_few_cells"
      return(list(
        status = status,
        all = data.frame(),
        significant = data.frame()
      ))
    }
    one <- tryCatch({
      obj <- subset(object, cells = cells)
      filtered <- .rc_drop_zero_count_atac_features(
        obj, atac_assay, paste0("Pando group ", group_id)
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
      infer_defaults <- list(object = grn, genes = target_genes)
      infer_defaults[names(pando_infer_args)] <- NULL
      grn <- do.call(Pando::infer_grn, c(infer_defaults, pando_infer_args))
      tab <- rc_extract_pando_tf_peak_gene(
        grn,
        sample_id = group_id,
        padj_threshold = padj_threshold,
        min_abs_estimate = min_abs_estimate,
        min_model_rsq = min_model_rsq,
        require_padj = require_padj
      )
      add_meta <- function(x) {
        if (!nrow(x)) return(x)
        x$group_id <- group_id
        for (column in group_cols) {
          x[[column]] <- as.character(vals[[column]][[1L]])
        }
        x[, c(
          "group_id", group_cols,
          setdiff(colnames(x), c("group_id", group_cols))
        ), drop = FALSE]
      }
      tab$all <- add_meta(tab$all)
      tab$significant <- add_meta(tab$significant)
      tab$peak_diagnostics <- filtered$diagnostics
      if (isTRUE(save_pando_objects)) {
        saveRDS(
          grn,
          file.path(
            outdir,
            "pando_objects",
            paste0(gsub("[^A-Za-z0-9_.-]+", "_", group_id), ".rds")
          )
        )
      }
      tab
    }, error = function(error) error)
    if (inherits(one, "error")) {
      status$status <- "failed"
      status$error_class <- class(one)[[1L]]
      status$error_message <- conditionMessage(one)
      if (identical(on_group_error, "stop")) stop(one)
      return(list(
        status = status,
        all = data.frame(),
        significant = data.frame()
      ))
    }
    status$n_atac_peaks_input <- one$peak_diagnostics$n_input_peaks
    status$n_zero_count_peaks_excluded <-
      one$peak_diagnostics$n_zero_count_peaks_excluded
    status$n_atac_peaks_used <- one$peak_diagnostics$n_retained_peaks
    status$status <- "ok"
    status$n_edges <- nrow(one$all)
    status$n_significant_edges <- nrow(one$significant)
    list(status = status, all = one$all, significant = one$significant)
  }

  results <- rc_parallel_lapply(group_ids, run_one_group, BPPARAM = BPPARAM)
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  all_edges <- do.call(rbind, lapply(results, `[[`, "all"))
  if (is.null(all_edges)) all_edges <- data.frame()
  significant <- do.call(rbind, lapply(results, `[[`, "significant"))
  if (is.null(significant)) significant <- data.frame()
  .rc_mm_write_tsv_gz(
    status, file.path(outdir, "pando_group_status.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    all_edges, file.path(outdir, "pando_tf_peak_gene_all.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    significant,
    file.path(outdir, "pando_tf_peak_gene_significant.tsv.gz")
  )
  failed <- status$status != "ok"
  if (any(failed)) {
    stop(
      paste0(
        "Every condition-by-cell-type Pando GRN must complete successfully. ",
        "Failed: ", paste(status$group_id[failed], collapse = "; ")
      ),
      call. = FALSE
    )
  }
  if (!nrow(significant)) {
    stop("No significant Pando TF-peak-gene edges were available.", call. = FALSE)
  }
  answer <- list(
    schema_version = "regcompass_single_cell_grn_v3",
    pando_installed_version = pando_install$version,
    pando_installation = pando_install,
    target_metabolic_genes = target_genes,
    sample_status = status,
    tf_peak_gene_all = all_edges,
    tf_peak_gene_significant = significant,
    normalization_policy = list(
      rna = "global single-cell NormalizeData before condition splitting",
      atac = paste(
        "cell-type-shared TF-IDF across conditions before condition",
        "splitting"
      ),
      zero_count_peaks = paste(
        "globally absent peaks are removed; cell-type-local absent peaks",
        "remain exact zeros and are not passed to RunTFIDF"
      ),
      pando_motifs = motif_policy,
      pando_regions = region_policy,
      pando_peak_cor = pando_infer_args$peak_cor %||% 0.01,
      pando_evidence_filters = list(
        min_cells = min_cells,
        padj_threshold = padj_threshold,
        min_abs_estimate = min_abs_estimate,
        min_model_rsq = min_model_rsq,
        require_padj = require_padj
      ),
      pando_rsq = paste0("finite rsq >= ", min_model_rsq)
    ),
    group_cols = group_cols
  )
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}
