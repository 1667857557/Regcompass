# Optional evidence-maximizing CORDA-like completion for Layer 2 union GEMs.

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
  model_completion <- match.arg(
    as.character(model_params$model_completion %||% "fastcore"),
    c("fastcore", "corda_like")
  )
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
    other_penalty = .rc_corda_scalar(
      model_params$corda_other_penalty %||% 1,
      "corda_other_penalty", 0, Inf
    ),
    negative_penalty = .rc_corda_scalar(
      model_params$corda_negative_penalty %||% 10,
      "corda_negative_penalty", 0, Inf
    ),
    include_evidence_outside_modules = .rc_corda_flag(
      model_params$corda_include_evidence_outside_modules %||% TRUE,
      "corda_include_evidence_outside_modules"
    ),
    max_medium_confidence_reactions = max_mc,
    evidence_definition = paste(
      "(1 - regulatory_weight) * max(within-cell-type RNA percentile,",
      "multiome percentile) + regulatory_weight * regulatory support"
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
  if (answer$negative_penalty < answer$other_penalty) {
    stop(
      "`corda_negative_penalty` must be greater than or equal to ",
      "`corda_other_penalty`.",
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
      "CORDA-like completion requires aligned Layer 1 matrices: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  reference <- layer1$reaction_expression
  for (name in required[-1L]) {
    if (!identical(dimnames(layer1[[name]]), dimnames(reference))) {
      stop("CORDA-like Layer 1 evidence matrices are not aligned.",
           call. = FALSE)
    }
  }
  params <- meta_modules$workflow_params
  celltype_col <- as.character(params$celltype_col %||% "cell_type")
  unit_meta <- layer1$unit_meta %||% layer1$metacell_meta
  if (!is.data.frame(unit_meta) || !celltype_col %in% colnames(unit_meta)) {
    stop("CORDA-like completion requires Layer 1 unit cell types.",
         call. = FALSE)
  }
  id_col <- if ("unit_id" %in% colnames(unit_meta)) {
    "unit_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    stop("CORDA-like completion requires unit_id or pool_id.", call. = FALSE)
  }
  unit_meta[[id_col]] <- as.character(unit_meta[[id_col]])
  unit_meta <- unit_meta[
    match(colnames(reference), unit_meta[[id_col]]), , drop = FALSE
  ]
  if (anyNA(unit_meta[[id_col]])) {
    stop("CORDA-like unit metadata do not align to Layer 1 columns.",
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
  answer$evidence_schema <- "regcompass_corda_like_reaction_evidence_v1"
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
      ordering <- order(
        -score[evidence_mc], evidence_mc, na.last = TRUE
      )
      evidence_mc <- evidence_mc[
        ordering[seq_len(max_medium_confidence_reactions)]
      ]
    }
  }
  hc <- core_reactions
  mc_module <- setdiff(module_reactions, hc)
  mc_evidence <- setdiff(evidence_mc, union(hc, mc_module))
  retained <- union(hc, union(mc_module, mc_evidence))
  rest <- setdiff(parent_reactions, retained)
  nc <- rest[
    is.finite(score[rest]) & score[rest] <= negative_confidence_threshold
  ]
  ot <- setdiff(rest, nc)
  class <- stats::setNames(rep("OT", length(parent_reactions)), parent_reactions)
  class[nc] <- "NC"
  class[mc_evidence] <- "MC_evidence"
  class[mc_module] <- "MC_module"
  class[hc] <- "HC"
  list(
    hc = hc,
    mc_module = mc_module,
    mc_evidence = mc_evidence,
    medium_confidence = union(mc_module, mc_evidence),
    biological = retained,
    ot = ot,
    nc = nc,
    evidence_score = score,
    evidence_class = class
  )
}

.rc_corda_support_costs <- function(
    reactions, classes, other_penalty = 1, negative_penalty = 10) {
  reactions <- unique(as.character(reactions))
  class <- as.character(classes$evidence_class[reactions])
  cost <- rep(other_penalty, length(reactions))
  cost[class == "NC"] <- negative_penalty
  cost[class %in% c("HC", "MC_module", "MC_evidence")] <- 0
  stats::setNames(cost, reactions)
}
