# RegCompass evidence adapter for original MATLAB CORDA2 HC/MC/NC/OT classes.

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
  if (length(value) != 1L || is.na(value) || !is.logical(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  isTRUE(value)
}

.rc_corda_rank_percentile <- function(value) {
  value <- suppressWarnings(as.numeric(value))
  answer <- rep(NA_real_, length(value))
  valid <- is.finite(value)
  if (!any(valid)) return(answer)
  ranked <- rank(value[valid], ties.method = "average", na.last = "keep")
  if (sum(valid) == 1L) {
    answer[valid] <- 1
  } else {
    answer[valid] <- (ranked - 1) / (sum(valid) - 1)
  }
  answer
}

.rc_corda_reaction_percentile <- function(value, reaction) {
  reaction <- as.character(reaction)
  value <- suppressWarnings(as.numeric(value))
  answer <- rep(NA_real_, length(value))
  groups <- split(seq_along(value), reaction)
  for (index in groups) answer[index] <- .rc_corda_rank_percentile(value[index])
  answer
}

.rc_corda_evidence_score <- function(
    rna_percentile, multiome_percentile, regulatory_support,
    regulatory_weight) {
  rna_percentile <- suppressWarnings(as.numeric(rna_percentile))
  multiome_percentile <- suppressWarnings(as.numeric(multiome_percentile))
  regulatory_support <- suppressWarnings(as.numeric(regulatory_support))
  regulatory_weight <- as.numeric(regulatory_weight)
  expression <- pmax(rna_percentile, multiome_percentile, na.rm = TRUE)
  expression[!is.finite(expression)] <- NA_real_
  answer <- (1 - regulatory_weight) * expression +
    regulatory_weight * regulatory_support
  answer[!is.finite(answer)] <- NA_real_
  pmin(1, pmax(0, answer))
}

.rc_corda_first_numeric_column <- function(tab, candidates) {
  candidate <- intersect(candidates, colnames(tab))
  if (!length(candidate)) return(rep(NA_real_, nrow(tab)))
  for (name in candidate) {
    value <- suppressWarnings(as.numeric(tab[[name]]))
    if (any(is.finite(value))) return(value)
  }
  rep(NA_real_, nrow(tab))
}

.rc_corda_layer1_support_table <- function(layer1) {
  candidates <- list(
    layer1$reaction_support,
    layer1$reaction_penalty,
    layer1$penalty_table,
    layer1$reaction_table,
    layer1$layer1_result$reaction_support,
    layer1$layer1_result$reaction_penalty,
    layer1$layer1_result$penalty_table,
    layer1$layer1_result$reaction_table
  )
  candidates <- candidates[vapply(candidates, is.data.frame, logical(1))]
  candidates <- candidates[vapply(candidates, function(tab) {
    all(c("reaction_id", "cell_type") %in% colnames(tab))
  }, logical(1))]
  if (!length(candidates)) {
    stop(
      "Layer 1 does not contain a reaction-support table with reaction_id ",
      "and cell_type columns required for CORDA2 evidence mapping.",
      call. = FALSE
    )
  }
  candidates[[1L]]
}

.rc_corda_layer1_regulatory_support <- function(layer1) {
  candidates <- list(
    layer1$condition_regulatory_support,
    layer1$regulatory_support,
    layer1$layer1_result$condition_regulatory_support,
    layer1$layer1_result$regulatory_support
  )
  candidates <- candidates[vapply(candidates, is.data.frame, logical(1))]
  candidates <- candidates[vapply(candidates, function(tab) {
    all(c("reaction_id", "cell_type") %in% colnames(tab))
  }, logical(1))]
  if (!length(candidates)) return(NULL)
  candidates[[1L]]
}

.rc_corda_layer1_multiome_support <- function(layer1) {
  candidates <- list(
    layer1$multiome_reaction_support,
    layer1$layer1_result$multiome_reaction_support,
    layer1$reaction_support,
    layer1$layer1_result$reaction_support
  )
  candidates <- candidates[vapply(candidates, is.data.frame, logical(1))]
  candidates <- candidates[vapply(candidates, function(tab) {
    all(c("reaction_id", "cell_type") %in% colnames(tab))
  }, logical(1))]
  if (!length(candidates)) return(NULL)
  candidates[[1L]]
}

.rc_corda_matrix_like <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

.rc_corda_native_layer1_root <- function(layer1) {
  candidates <- list(layer1, layer1$layer1_result)
  candidates <- candidates[vapply(candidates, is.list, logical(1))]
  for (candidate in candidates) {
    if (.rc_corda_matrix_like(candidate$reaction_expression_rna_only) &&
        .rc_corda_matrix_like(candidate$reaction_expression) &&
        is.data.frame(candidate$unit_meta %||% candidate$metacell_meta)) {
      return(candidate)
    }
  }
  NULL
}

.rc_corda_row_medians <- function(x) {
  value <- matrixStats::rowMedians(as.matrix(x), na.rm = TRUE)
  value[!is.finite(value)] <- NA_real_
  value
}

.rc_corda_row_maxs <- function(x) {
  value <- matrixStats::rowMaxs(as.matrix(x), na.rm = TRUE)
  value[!is.finite(value)] <- NA_real_
  value
}

.rc_corda_condition_median_max <- function(
    value, unit_meta, condition, cell_type) {
  reactions <- rownames(value)
  cell_types <- sort(unique(cell_type))
  rows <- lapply(cell_types, function(current_cell_type) {
    selected <- which(cell_type == current_cell_type)
    current_conditions <- sort(unique(condition[selected]))
    condition_medians <- lapply(current_conditions, function(current_condition) {
      index <- which(
        cell_type == current_cell_type & condition == current_condition
      )
      .rc_corda_row_medians(value[, index, drop = FALSE])
    })
    condition_matrix <- do.call(cbind, condition_medians)
    data.frame(
      cell_type = current_cell_type,
      reaction_id = reactions,
      value = .rc_corda_row_maxs(condition_matrix),
      n_conditions = length(current_conditions),
      n_metacells = length(selected),
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer
}

.rc_corda_native_layer1_evidence_tables <- function(layer1) {
  root <- .rc_corda_native_layer1_root(layer1)
  if (is.null(root)) return(NULL)

  rna <- root$reaction_expression_rna_only
  multiome <- root$reaction_expression
  regulatory <- root$reaction_regulatory_support_fraction
  unit_meta <- as.data.frame(root$unit_meta %||% root$metacell_meta)
  params <- root$workflow_params %||% layer1$workflow_params %||% list()

  if (!identical(dimnames(rna), dimnames(multiome))) {
    stop(
      "Layer 1 RNA-only and multiome reaction matrices must align exactly ",
      "for CORDA2 evidence mapping.",
      call. = FALSE
    )
  }
  if (is.null(rownames(rna)) || is.null(colnames(rna)) ||
      anyNA(rownames(rna)) || any(!nzchar(rownames(rna))) ||
      anyDuplicated(rownames(rna)) || anyNA(colnames(rna)) ||
      any(!nzchar(colnames(rna))) || anyDuplicated(colnames(rna))) {
    stop(
      "Layer 1 reaction matrices require unique non-empty reaction and ",
      "metacell identifiers for CORDA2 evidence mapping.",
      call. = FALSE
    )
  }
  if (!is.null(regulatory)) {
    if (!.rc_corda_matrix_like(regulatory) ||
        !identical(dimnames(regulatory), dimnames(rna))) {
      stop(
        "Layer 1 reaction regulatory support must align exactly with the ",
        "reaction matrices.",
        call. = FALSE
      )
    }
  }

  unit_col <- intersect(
    c("unit_id", "metacell_id", "pool_id"), colnames(unit_meta)
  )
  if (!length(unit_col)) {
    stop(
      "Layer 1 unit metadata require unit_id, metacell_id or pool_id for ",
      "CORDA2 evidence mapping.",
      call. = FALSE
    )
  }
  unit_col <- unit_col[[1L]]
  unit_id <- as.character(unit_meta[[unit_col]])
  if (anyNA(unit_id) || any(!nzchar(unit_id)) || anyDuplicated(unit_id)) {
    stop("Layer 1 unit identifiers must be unique and non-empty.",
         call. = FALSE)
  }
  index <- match(colnames(rna), unit_id)
  if (anyNA(index)) {
    stop(
      "Layer 1 reaction matrix columns are absent from unit metadata.",
      call. = FALSE
    )
  }
  unit_meta <- unit_meta[index, , drop = FALSE]

  celltype_col <- as.character(params$celltype_col %||% "cell_type")
  if (length(celltype_col) != 1L || is.na(celltype_col) ||
      !nzchar(celltype_col) || !celltype_col %in% colnames(unit_meta)) {
    stop(
      "Layer 1 workflow parameters do not resolve a valid cell-type column ",
      "for CORDA2 evidence mapping.",
      call. = FALSE
    )
  }
  condition_col <- params$condition_col
  if (is.null(condition_col) || !length(condition_col) ||
      is.na(condition_col[[1L]]) || !nzchar(as.character(condition_col[[1L]]))) {
    condition <- rep("all", nrow(unit_meta))
    condition_col <- NA_character_
  } else {
    condition_col <- as.character(condition_col[[1L]])
    if (!condition_col %in% colnames(unit_meta)) {
      stop(
        "Layer 1 workflow parameters name a condition column absent from ",
        "unit metadata.",
        call. = FALSE
      )
    }
    condition <- trimws(as.character(unit_meta[[condition_col]]))
  }
  cell_type <- trimws(as.character(unit_meta[[celltype_col]]))
  if (anyNA(cell_type) || any(!nzchar(cell_type)) ||
      anyNA(condition) || any(!nzchar(condition))) {
    stop(
      "Layer 1 metacells require complete cell-type and condition labels for ",
      "CORDA2 evidence mapping.",
      call. = FALSE
    )
  }

  rna_table <- .rc_corda_condition_median_max(
    rna, unit_meta, condition, cell_type
  )
  names(rna_table)[names(rna_table) == "value"] <- "rna_reaction_support"
  multiome_table <- .rc_corda_condition_median_max(
    multiome, unit_meta, condition, cell_type
  )
  names(multiome_table)[names(multiome_table) == "value"] <-
    "multiome_reaction_support"

  regulatory_table <- NULL
  if (!is.null(regulatory)) {
    regulatory_table <- .rc_corda_condition_median_max(
      regulatory, unit_meta, condition, cell_type
    )
    names(regulatory_table)[names(regulatory_table) == "value"] <-
      "regulatory_support"
  }

  list(
    support = rna_table,
    multiome = multiome_table,
    regulatory = regulatory_table,
    source = "native_layer1_metacell_matrices",
    aggregation = "condition_median_max",
    aggregation_contract = paste(
      "median across metacells within each condition and cell type, then max",
      "across conditions within the same cell type"
    ),
    condition_col = condition_col,
    celltype_col = celltype_col
  )
}

.rc_corda_layer1_evidence_tables <- function(layer1) {
  native <- .rc_corda_native_layer1_evidence_tables(layer1)
  if (!is.null(native)) return(native)
  list(
    support = .rc_corda_layer1_support_table(layer1),
    multiome = .rc_corda_layer1_multiome_support(layer1),
    regulatory = .rc_corda_layer1_regulatory_support(layer1),
    source = "legacy_reaction_support_tables",
    aggregation = "legacy_celltype_mean",
    aggregation_contract =
      "legacy reaction-support rows aggregated by cell-type mean"
  )
}

.rc_corda_aggregate_support <- function(
    tab, value, value_name, aggregation = c("max", "mean")) {
  aggregation <- match.arg(aggregation)
  value <- suppressWarnings(as.numeric(value))
  key <- interaction(
    as.character(tab$cell_type), as.character(tab$reaction_id),
    drop = TRUE, lex.order = TRUE
  )
  rows <- lapply(split(seq_len(nrow(tab)), key), function(index) {
    values <- value[index]
    values <- values[is.finite(values)]
    aggregate_value <- if (!length(values)) {
      NA_real_
    } else if (identical(aggregation, "max")) {
      max(values)
    } else {
      mean(values)
    }
    data.frame(
      cell_type = as.character(tab$cell_type[index[[1L]]]),
      reaction_id = as.character(tab$reaction_id[index[[1L]]]),
      value = aggregate_value,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  names(answer)[names(answer) == "value"] <- value_name
  rownames(answer) <- NULL
  answer
}

.rc_corda_reaction_evidence <- function(
    layer1, meta_modules, regulatory_weight = 0.20) {
  tables <- .rc_corda_layer1_evidence_tables(layer1)
  support <- tables$support
  support$cell_type <- as.character(support$cell_type)
  support$reaction_id <- as.character(support$reaction_id)
  support <- support[
    !is.na(support$cell_type) & nzchar(support$cell_type) &
      !is.na(support$reaction_id) & nzchar(support$reaction_id),
    , drop = FALSE
  ]
  if (!nrow(support)) {
    stop("Layer 1 reaction support is empty after validation.", call. = FALSE)
  }

  rna <- .rc_corda_first_numeric_column(
    support,
    c(
      "rna_reaction_support", "rna_support", "expression_support",
      "reaction_support", "support", "capacity", "penalty"
    )
  )
  if ("penalty" %in% colnames(support) &&
      all(!is.finite(rna) | rna == suppressWarnings(as.numeric(support$penalty)))) {
    penalty <- suppressWarnings(as.numeric(support$penalty))
    rna <- ifelse(is.finite(penalty), -penalty, NA_real_)
  }
  rna_table <- .rc_corda_aggregate_support(
    support, rna, "rna_support", aggregation = "mean"
  )
  rna_table$rna_percentile <- .rc_corda_reaction_percentile(
    rna_table$rna_support, rna_table$cell_type
  )

  multiome_tab <- tables$multiome
  if (is.null(multiome_tab)) {
    multiome_table <- rna_table[
      c("cell_type", "reaction_id", "rna_support", "rna_percentile")
    ]
    names(multiome_table)[3:4] <- c(
      "multiome_support", "multiome_percentile"
    )
  } else {
    multiome_tab$cell_type <- as.character(multiome_tab$cell_type)
    multiome_tab$reaction_id <- as.character(multiome_tab$reaction_id)
    multiome <- .rc_corda_first_numeric_column(
      multiome_tab,
      c(
        "multiome_reaction_support", "multiome_support",
        "regulatory_expression_support", "reaction_support", "support"
      )
    )
    multiome_table <- .rc_corda_aggregate_support(
      multiome_tab, multiome, "multiome_support", aggregation = "mean"
    )
    multiome_table$multiome_percentile <- .rc_corda_reaction_percentile(
      multiome_table$multiome_support, multiome_table$cell_type
    )
  }

  regulatory_tab <- tables$regulatory
  if (is.null(regulatory_tab)) {
    regulatory_table <- rna_table[c("cell_type", "reaction_id")]
    regulatory_table$regulatory_support <- 0
  } else {
    regulatory_tab$cell_type <- as.character(regulatory_tab$cell_type)
    regulatory_tab$reaction_id <- as.character(regulatory_tab$reaction_id)
    regulatory <- .rc_corda_first_numeric_column(
      regulatory_tab,
      c(
        "regulatory_support", "condition_regulatory_support",
        "tf_atac_support", "edge_support", "support"
      )
    )
    regulatory_table <- .rc_corda_aggregate_support(
      regulatory_tab, regulatory, "regulatory_support", aggregation = "max"
    )
    regulatory_table$regulatory_support <- pmin(
      1, pmax(0, regulatory_table$regulatory_support)
    )
    regulatory_table$regulatory_support[
      !is.finite(regulatory_table$regulatory_support)
    ] <- 0
  }

  answer <- merge(
    rna_table, multiome_table,
    by = c("cell_type", "reaction_id"), all = TRUE, sort = FALSE
  )
  answer <- merge(
    answer, regulatory_table,
    by = c("cell_type", "reaction_id"), all = TRUE, sort = FALSE
  )
  answer$regulatory_support[
    !is.finite(answer$regulatory_support)
  ] <- 0
  answer$evidence_score <- .rc_corda_evidence_score(
    answer$rna_percentile,
    answer$multiome_percentile,
    answer$regulatory_support,
    regulatory_weight = regulatory_weight
  )

  membership <- .rc_meta_module_reaction_membership(meta_modules)
  if (!"cell_type" %in% colnames(membership)) {
    stop("Meta-module reaction membership lacks cell_type.", call. = FALSE)
  }
  membership$cell_type <- as.character(membership$cell_type)
  membership$reaction_id <- as.character(membership$reaction_id)
  module_flag <- unique(membership[c("cell_type", "reaction_id")])
  module_flag$merged_meta_module_member <- TRUE
  answer <- merge(
    answer, module_flag,
    by = c("cell_type", "reaction_id"), all = TRUE, sort = FALSE
  )
  answer$merged_meta_module_member[
    is.na(answer$merged_meta_module_member)
  ] <- FALSE
  answer$evidence_source <- tables$source
  answer$evidence_aggregation <- tables$aggregation
  answer$evidence_aggregation_contract <- tables$aggregation_contract
  answer$evidence_contract <- paste(
    "RegCompass maps multiome evidence to original CORDA2 HC, MC, NC and OT",
    "reaction groups before the reconstruction state machine begins; native",
    "Layer 1 metacells use condition_median_max aggregation"
  )
  answer
}

.rc_layer2_corda_reaction_evidence <- function(
    layer1, meta_modules, regulatory_weight = 0.20) {
  .rc_corda_reaction_evidence(
    layer1 = layer1,
    meta_modules = meta_modules,
    regulatory_weight = regulatory_weight
  )
}

.rc_corda_classify_reactions <- function(
    parent_reactions, module_reactions, core_reactions,
    reaction_evidence, medium_confidence_threshold,
    negative_confidence_threshold,
    include_evidence_outside_modules,
    max_medium_confidence_reactions = Inf) {
  reactions <- unique(as.character(parent_reactions))
  evidence <- reaction_evidence[
    match(reactions, as.character(reaction_evidence$reaction_id)),
    , drop = FALSE
  ]
  score <- suppressWarnings(as.numeric(evidence$evidence_score))
  names(score) <- reactions
  module_reactions <- intersect(
    unique(as.character(module_reactions)), reactions
  )
  core_reactions <- intersect(
    unique(as.character(core_reactions)), module_reactions
  )
  mc_module <- setdiff(module_reactions, core_reactions)
  outside <- setdiff(reactions, module_reactions)
  mc_evidence <- character()
  if (isTRUE(include_evidence_outside_modules)) {
    mc_evidence <- outside[
      is.finite(score[outside]) &
        score[outside] >= medium_confidence_threshold
    ]
    if (is.finite(max_medium_confidence_reactions) &&
        length(mc_evidence) > max_medium_confidence_reactions) {
      order_index <- order(
        -score[mc_evidence], mc_evidence, method = "radix"
      )
      mc_evidence <- mc_evidence[
        order_index[seq_len(max_medium_confidence_reactions)]
      ]
    }
  }
  mc <- union(mc_module, mc_evidence)
  remaining <- setdiff(reactions, union(core_reactions, mc))
  nc <- remaining[
    is.finite(score[remaining]) &
      score[remaining] <= negative_confidence_threshold
  ]
  ot <- setdiff(remaining, nc)

  confidence <- stats::setNames(rep("OT", length(reactions)), reactions)
  confidence[core_reactions] <- "HC"
  confidence[mc_module] <- "MC_module"
  confidence[mc_evidence] <- "MC_evidence"
  confidence[nc] <- "NC"
  list(
    confidence = confidence,
    initial_confidence = confidence,
    hc = core_reactions,
    mc_module = mc_module,
    mc_evidence = mc_evidence,
    mc = mc,
    nc = nc,
    ot = ot,
    evidence_score = score,
    confidence_contract = paste(
      "HC=merged core; MC=remaining merged module plus selected high-evidence",
      "outside-module reactions; NC=finite low-evidence; OT=remaining"
    )
  )
}
