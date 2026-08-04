# Cell-type reaction confidence mapping for the original CORDA algorithm.

.rc_corda_scalar <- function(value, name, lower = -Inf, upper = Inf,
                             finite = TRUE) {
  value <- suppressWarnings(as.numeric(value))
  valid <- length(value) == 1L && !is.na(value) &&
    (!finite || is.finite(value)) && value >= lower && value <= upper
  if (!valid) {
    stop("`", name, "` must be one number in [", lower, ", ", upper,
         "].", call. = FALSE)
  }
  value
}

.rc_corda_integer <- function(value, name, lower = 0L) {
  numeric_value <- suppressWarnings(as.numeric(value))
  integer_value <- suppressWarnings(as.integer(numeric_value))
  if (length(numeric_value) != 1L || is.na(numeric_value) ||
      !is.finite(numeric_value) || numeric_value != integer_value ||
      integer_value < lower) {
    stop("`", name, "` must be one integer >= ", lower, ".",
         call. = FALSE)
  }
  integer_value
}

.rc_corda_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(value)
}

.rc_layer2_corda_options <- function(model_params = list()) {
  if (!is.list(model_params)) {
    stop("`layer2_args$model_params` must be a list.", call. = FALSE)
  }
  requested <- as.character(model_params$model_completion %||% "fastcore")
  if (length(requested) != 1L || is.na(requested)) {
    stop("`model_completion` must be `fastcore` or `corda`.", call. = FALSE)
  }
  if (identical(requested, "corda_like")) requested <- "corda"
  model_completion <- match.arg(requested, c("fastcore", "corda"))
  obsolete <- intersect(
    names(model_params),
    c("corda_other_penalty", "corda_negative_penalty")
  )
  if (length(obsolete)) {
    stop(
      "Obsolete weighted-FASTCORE parameter(s): ",
      paste(obsolete, collapse = ", "),
      ". Original CORDA uses `corda_gamma`, `corda_kappa`, `corda_n`, ",
      "and `corda_p`.",
      call. = FALSE
    )
  }
  max_mc <- suppressWarnings(as.numeric(
    model_params$corda_max_medium_confidence_reactions %||% Inf
  ))
  if (length(max_mc) != 1L || is.na(max_mc) || max_mc < 0) {
    stop(
      "`corda_max_medium_confidence_reactions` must be non-negative or Inf.",
      call. = FALSE
    )
  }
  if (is.finite(max_mc)) max_mc <- as.integer(floor(max_mc))
  answer <- list(
    model_completion = model_completion,
    gamma = .rc_corda_scalar(
      model_params$corda_gamma %||% 1e5,
      "corda_gamma", 1, Inf
    ),
    kappa = .rc_corda_scalar(
      model_params$corda_kappa %||% 1e-2,
      "corda_kappa", 0, Inf
    ),
    epsilon = .rc_corda_scalar(
      model_params$corda_epsilon %||% 1,
      "corda_epsilon", .Machine$double.eps, Inf
    ),
    n = .rc_corda_integer(
      model_params$corda_n %||% 5L,
      "corda_n", 1L
    ),
    p = .rc_corda_integer(
      model_params$corda_p %||% 2L,
      "corda_p", 1L
    ),
    seed = .rc_corda_integer(
      model_params$corda_seed %||% 1L,
      "corda_seed", 0L
    ),
    flux_tolerance = .rc_corda_scalar(
      model_params$corda_flux_tolerance %||% 1e-8,
      "corda_flux_tolerance", .Machine$double.eps, Inf
    ),
    medium_confidence_threshold = .rc_corda_scalar(
      model_params$corda_medium_confidence_threshold %||% 0.75,
      "corda_medium_confidence_threshold", 0, 1
    ),
    negative_confidence_threshold = .rc_corda_scalar(
      model_params$corda_negative_confidence_threshold %||% 0.10,
      "corda_negative_confidence_threshold", 0, 1
    ),
    regulatory_weight = .rc_corda_scalar(
      model_params$corda_regulatory_weight %||% 0.20,
      "corda_regulatory_weight", 0, 1
    ),
    include_evidence_outside_modules = .rc_corda_flag(
      model_params$corda_include_evidence_outside_modules %||% TRUE,
      "corda_include_evidence_outside_modules"
    ),
    max_medium_confidence_reactions = max_mc,
    evidence_definition = paste(
      "(1 - regulatory_weight) * max(within-cell-type RNA percentile,",
      "multiome percentile) + regulatory_weight * regulatory support"
    ),
    algorithm = "Schultz_Qutub_CORDA_2016_three_stage_dependency_assessment",
    paper_defaults = c(
      gamma = 1e5, kappa = 1e-2, epsilon = 1, n = 5, p = 2
    )
  )
  if (answer$negative_confidence_threshold >
      answer$medium_confidence_threshold) {
    stop(
      "`corda_negative_confidence_threshold` must not exceed ",
      "`corda_medium_confidence_threshold`.",
      call. = FALSE
    )
  }
  answer
}

.rc_corda_rank01 <- function(value) {
  value <- as.numeric(value)
  answer <- rep(NA_real_, length(value))
  keep <- is.finite(value)
  n <- sum(keep)
  if (!n) return(answer)
  if (n == 1L) {
    answer[keep] <- 1
  } else {
    answer[keep] <- (rank(value[keep], ties.method = "average") - 1) /
      (n - 1)
  }
  answer
}

.rc_corda_row_summary <- function(x, units) {
  x <- as.matrix(x[, units, drop = FALSE])
  apply(x, 1L, function(value) {
    value <- value[is.finite(value)]
    if (length(value)) stats::median(value) else NA_real_
  })
}

.rc_layer2_corda_reaction_evidence <- function(
    layer1, meta_modules, regulatory_weight = 0.20) {
  required <- c(
    "reaction_expression", "reaction_expression_rna_only",
    "reaction_regulatory_support_fraction"
  )
  missing <- required[!vapply(required, function(name) {
    value <- layer1[[name]]
    is.numeric(value) && !is.null(dim(value))
  }, logical(1))]
  if (length(missing)) {
    stop(
      "CORDA completion requires aligned Layer 1 matrices: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  reference <- layer1$reaction_expression
  for (name in required[-1L]) {
    if (!identical(dimnames(layer1[[name]]), dimnames(reference))) {
      stop("CORDA Layer 1 evidence matrices are not aligned.",
           call. = FALSE)
    }
  }
  params <- meta_modules$workflow_params
  celltype_col <- as.character(params$celltype_col %||% "cell_type")
  unit_meta <- layer1$unit_meta %||% layer1$metacell_meta
  if (!is.data.frame(unit_meta) || !celltype_col %in% colnames(unit_meta)) {
    stop("CORDA completion requires Layer 1 unit cell types.",
         call. = FALSE)
  }
  id_col <- if ("unit_id" %in% colnames(unit_meta)) {
    "unit_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("CORDA completion requires unit_id or pool_id.", call. = FALSE)
  }
  unit_meta[[id_col]] <- as.character(unit_meta[[id_col]])
  unit_meta <- unit_meta[
    match(colnames(reference), unit_meta[[id_col]]), , drop = FALSE
  ]
  if (anyNA(unit_meta[[id_col]])) {
    stop("CORDA unit metadata do not align to Layer 1 columns.",
         call. = FALSE)
  }
  rows <- lapply(unique(as.character(unit_meta[[celltype_col]])), function(ct) {
    units <- unit_meta[[id_col]][as.character(unit_meta[[celltype_col]]) == ct]
    rna <- .rc_corda_row_summary(layer1$reaction_expression_rna_only, units)
    multiome <- .rc_corda_row_summary(layer1$reaction_expression, units)
    regulatory <- .rc_corda_row_summary(
      layer1$reaction_regulatory_support_fraction, units
    )
    rna_percentile <- .rc_corda_rank01(rna)
    multiome_percentile <- .rc_corda_rank01(multiome)
    expression_percentile <- pmax(
      rna_percentile, multiome_percentile, na.rm = TRUE
    )
    expression_percentile[
      !is.finite(rna_percentile) & !is.finite(multiome_percentile)
    ] <- NA_real_
    regulatory_clamped <- pmin(pmax(regulatory, 0), 1)
    regulatory_clamped[!is.finite(regulatory_clamped)] <- 0
    evidence_score <-
      (1 - regulatory_weight) * expression_percentile +
      regulatory_weight * regulatory_clamped
    evidence_score[!is.finite(expression_percentile)] <- NA_real_
    data.frame(
      cell_type = ct,
      reaction_id = rownames(reference),
      rna_median = as.numeric(rna),
      multiome_median = as.numeric(multiome),
      rna_percentile = as.numeric(rna_percentile),
      multiome_percentile = as.numeric(multiome_percentile),
      regulatory_support = as.numeric(regulatory),
      evidence_score = as.numeric(evidence_score),
      stringsAsFactors = FALSE
    )
  })
  answer <- .rc_bind_frames_fill(rows)
  answer$evidence_schema <- "regcompass_corda_reaction_evidence_v2"
  answer
}

.rc_corda_classify_reactions <- function(
    parent_reactions, module_reactions, core_reactions, reaction_evidence,
    medium_confidence_threshold = 0.75,
    negative_confidence_threshold = 0.10,
    include_evidence_outside_modules = TRUE,
    max_medium_confidence_reactions = Inf) {
  parent_reactions <- unique(as.character(parent_reactions))
  module_reactions <- intersect(
    unique(as.character(module_reactions)), parent_reactions
  )
  core_reactions <- intersect(
    unique(as.character(core_reactions)), module_reactions
  )
  evidence <- reaction_evidence[
    match(parent_reactions, as.character(reaction_evidence$reaction_id)),
    , drop = FALSE
  ]
  evidence$reaction_id <- parent_reactions
  score <- suppressWarnings(as.numeric(evidence$evidence_score))
  names(score) <- parent_reactions
  outside <- setdiff(parent_reactions, module_reactions)
  evidence_mc <- character()
  if (isTRUE(include_evidence_outside_modules)) {
    evidence_mc <- outside[
      is.finite(score[outside]) &
        score[outside] >= medium_confidence_threshold
    ]
    if (is.finite(max_medium_confidence_reactions) &&
        length(evidence_mc) > max_medium_confidence_reactions) {
      ordering <- order(-score[evidence_mc], evidence_mc, na.last = TRUE)
      evidence_mc <- evidence_mc[
        ordering[seq_len(max_medium_confidence_reactions)]
      ]
    }
  }
  hc <- core_reactions
  mc_module <- setdiff(module_reactions, hc)
  mc_evidence <- setdiff(evidence_mc, union(hc, mc_module))
  mc <- union(mc_module, mc_evidence)
  remaining <- setdiff(parent_reactions, union(hc, mc))
  nc <- remaining[
    is.finite(score[remaining]) &
      score[remaining] <= negative_confidence_threshold
  ]
  ot <- setdiff(remaining, nc)
  confidence <- stats::setNames(rep("OT", length(parent_reactions)),
                                parent_reactions)
  confidence[nc] <- "NC"
  confidence[mc_evidence] <- "MC_evidence"
  confidence[mc_module] <- "MC_module"
  confidence[hc] <- "HC"
  list(
    hc = hc,
    mc_module = mc_module,
    mc_evidence = mc_evidence,
    mc = mc,
    nc = nc,
    ot = ot,
    evidence_score = score,
    confidence = confidence,
    initial_confidence = confidence,
    confidence_contract = paste(
      "HC=merged core; MC=non-core module plus optional high-evidence",
      "outside-module reactions; NC=finite low evidence; OT=remaining"
    )
  )
}
