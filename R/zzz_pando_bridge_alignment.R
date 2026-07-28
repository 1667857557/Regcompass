# Final Pando-to-RegCompass bridge validation and execution semantics.

.rc_validate_pando_bridge_list <- function(x, label) {
  if (!is.list(x)) {
    stop("`", label, "` must be a list.", call. = FALSE)
  }
  invisible(TRUE)
}

.rc_reject_reserved_pando_args <- function(x, reserved, label) {
  conflict <- intersect(names(x), reserved)
  if (length(conflict)) {
    stop(
      "`", label, "` cannot override RegCompass-managed fields: ",
      paste(conflict, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.rc_validate_pando_bridge_args <- function(
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list()) {
  .rc_validate_pando_bridge_list(pando_initiate_args, "pando_initiate_args")
  .rc_validate_pando_bridge_list(pando_motif_args, "pando_motif_args")
  .rc_validate_pando_bridge_list(pando_infer_args, "pando_infer_args")

  .rc_reject_reserved_pando_args(
    pando_initiate_args,
    c("object", "peak_assay", "rna_assay"),
    "pando_initiate_args"
  )
  .rc_reject_reserved_pando_args(
    pando_motif_args,
    c("object", "pfm", "genome"),
    "pando_motif_args"
  )
  .rc_reject_reserved_pando_args(
    pando_infer_args,
    c(
      "object", "cell_type_col", "condition_col", "genes", "network_name",
      "min_cells_per_condition", "on_small_condition", "BPPARAM"
    ),
    "pando_infer_args"
  )

  aggregate_fields <- intersect(
    names(pando_infer_args),
    c("aggregate_rna_col", "aggregate_peaks_col")
  )
  if (length(aggregate_fields)) {
    stop(
      "Canonical RegCompass Stage 1 fits Pando on paired single cells. ",
      "Do not supply ", paste(aggregate_fields, collapse = ", "),
      "; RegCompass Stage 2 owns metacell aggregation.",
      call. = FALSE
    )
  }

  if (!is.null(pando_infer_args$candidate_screen)) {
    candidate_screen <- as.character(pando_infer_args$candidate_screen)
    if (length(candidate_screen) != 1L || is.na(candidate_screen) ||
        !candidate_screen %in% c("motif_domain", "condition_union", "pooled")) {
      stop(
        "`pando_infer_args$candidate_screen` must be one of ",
        "\"motif_domain\", \"condition_union\", or \"pooled\".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.rc_validate_pando_repository_before_condition_api <-
  .rc_validate_pando_repository

.rc_validate_pando_repository <- function(
    description = NULL, installed_version = NULL,
    validate_api = is.null(description)) {
  result <- .rc_validate_pando_repository_before_condition_api(
    description = description,
    installed_version = installed_version,
    validate_api = validate_api
  )
  result$condition_api_verified <- FALSE
  if (isTRUE(validate_api)) {
    namespace <- asNamespace("Pando")
    required_condition_api <- c("infer_condition_grn", "condition_grn_fit")
    missing <- required_condition_api[!vapply(
      required_condition_api,
      function(name) exists(name, envir = namespace, inherits = FALSE),
      logical(1)
    )]
    if (length(missing)) {
      stop(
        "Installed Pando is missing the condition-comparable API: ",
        paste(missing, collapse = ", "),
        ". Install 1667857557/Pando_regcompass >= 1.2.1.",
        call. = FALSE
      )
    }
    if (utils::package_version(result$version) <
        utils::package_version("1.2.1")) {
      stop(
        "RegCompass requires Pando >= 1.2.1 for the explicit ",
        "ConditionGRNFit comparison mask; installed version is ",
        result$version, ".",
        call. = FALSE
      )
    }
    result$condition_api_verified <- TRUE
  }
  result
}

.rc_condition_fit_comparison_mask <- function(fit) {
  beta <- as.matrix(fit$beta)
  eligibility <- as.matrix(fit$eligibility_mask)
  if (!is.logical(eligibility) || anyNA(eligibility) ||
      !identical(dim(eligibility), dim(beta)) ||
      !identical(dimnames(eligibility), dimnames(beta))) {
    stop("ConditionGRNFit eligibility mask is invalid.", call. = FALSE)
  }
  if (is.null(fit$comparison_mask)) {
    stop(
      "ConditionGRNFit lacks the explicit `comparison_mask`; install ",
      "1667857557/Pando_regcompass >= 1.2.1.",
      call. = FALSE
    )
  }
  comparison <- as.matrix(fit$comparison_mask)
  if (!is.logical(comparison) || anyNA(comparison) ||
      !identical(dim(comparison), dim(beta)) ||
      !identical(dimnames(comparison), dimnames(beta))) {
    stop("ConditionGRNFit comparison mask is invalid.", call. = FALSE)
  }
  if (!fit$reference_condition %in% colnames(eligibility)) {
    stop("ConditionGRNFit reference condition is absent.", call. = FALSE)
  }
  expected <- eligibility & matrix(
    eligibility[, fit$reference_condition],
    nrow = nrow(eligibility),
    ncol = ncol(eligibility)
  )
  dimnames(expected) <- dimnames(eligibility)
  if (!identical(comparison, expected)) {
    stop(
      "ConditionGRNFit comparison mask is inconsistent with eligibility.",
      call. = FALSE
    )
  }
  comparison
}

.rc_run_condition_single_cell_grns_before_bridge_alignment <-
  .rc_run_condition_single_cell_grns

.rc_run_condition_single_cell_grns <- function(
    object, gem, outdir, pfm = NULL, genome,
    condition_col = "condition", celltype_col = "cell_type",
    rna_assay = "RNA", atac_assay = "ATAC",
    min_cells = 20L,
    pando_initiate_args = list(exclude_exons = TRUE),
    pando_motif_args = list(),
    pando_infer_args = list(
      method = "shared_design_independent",
      candidate_screen = "motif_domain",
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
  .rc_validate_pando_bridge_args(
    pando_initiate_args = pando_initiate_args,
    pando_motif_args = pando_motif_args,
    pando_infer_args = pando_infer_args
  )
  if (is.null(pando_infer_args$candidate_screen)) {
    pando_infer_args$candidate_screen <- "motif_domain"
  }
  result <- .rc_run_condition_single_cell_grns_before_bridge_alignment(
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
    pando_infer_args = pando_infer_args,
    padj_threshold = padj_threshold,
    min_abs_estimate = min_abs_estimate,
    min_model_rsq = min_model_rsq,
    require_padj = require_padj,
    save_pando_objects = save_pando_objects,
    BPPARAM = BPPARAM,
    on_group_error = on_group_error,
    species = species
  )
  result$comparison_policy <- paste(
    "condition effects require Pando comparison_mask support in both",
    "the requested condition and the explicit reference condition"
  )
  result$normalization_policy$pando_candidate_screen <-
    pando_infer_args$candidate_screen
  result$normalization_policy$comparison_support <- result$comparison_policy
  saveRDS(result, file.path(outdir, "single_cell_grn.rds"))
  result
}

.rc_regcompass_step_grn_before_parallel_alignment <-
  rc_regcompass_step_grn

rc_regcompass_step_grn <- function(
    object, gem, outdir, genome,
    pfm = NULL,
    species = c("auto", "human", "mouse"),
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    pando_args = list(),
    parallel = TRUE,
    BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  if (!is.list(pando_args)) {
    stop("`pando_args` must be a list.", call. = FALSE)
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  if (identical(BPPARAM, TRUE)) {
    stop(
      "`BPPARAM = TRUE` is invalid. Supply a BiocParallelParam object, ",
      "`NULL`, or `FALSE`.",
      call. = FALSE
    )
  }
  supplied_bpparam <- !is.null(BPPARAM) && !identical(BPPARAM, FALSE)
  if (!isTRUE(parallel) && supplied_bpparam) {
    warning(
      "Ignoring `BPPARAM` because `parallel = FALSE`.",
      call. = FALSE
    )
  }
  use_bpparam <- isTRUE(parallel) && supplied_bpparam
  use_pando_native <- isTRUE(parallel) && !use_bpparam

  infer_args <- pando_args$pando_infer_args %||% list()
  if (!is.list(infer_args)) {
    stop("`pando_args$pando_infer_args` must be a list.", call. = FALSE)
  }
  if (!is.null(infer_args$parallel) &&
      !identical(isTRUE(infer_args$parallel), use_pando_native)) {
    warning(
      "Stage-level `parallel` and `BPPARAM` control Pando execution; ",
      "ignoring `pando_args$pando_infer_args$parallel`.",
      call. = FALSE
    )
  }
  infer_args$parallel <- use_pando_native
  pando_args$pando_infer_args <- infer_args

  result <- .rc_regcompass_step_grn_before_parallel_alignment(
    object = object,
    gem = gem,
    outdir = outdir,
    genome = genome,
    pfm = pfm,
    species = species,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    pando_args = pando_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )
  result$params$pando_parallel <- list(
    requested = parallel,
    backend = if (use_bpparam) {
      paste(class(BPPARAM), collapse = "/")
    } else if (use_pando_native) {
      "Pando_native"
    } else {
      "serial"
    },
    infer_parallel = use_pando_native,
    BPPARAM_supplied = use_bpparam
  )
  saveRDS(result, file.path(outdir, "step_grn.rds"))
  result
}
