# Workflow stage 1: establish the canonical biological and inference contracts.
# This file is explicitly collated after the base implementations so existing
# APIs remain available while the integrated workflow uses safer defaults.

.rc_original_reaction_capacity <- rc_reaction_capacity
.rc_original_prepare_humangem_gpr_table <- rc_prepare_humangem_gpr_table
.rc_original_project_metabolic_grn <- rc_project_metabolic_grn
.rc_original_compute_multiome_penalty_core <- .rc_compute_multiome_penalty_core
.rc_original_layer2_penalty <- rc_layer2_penalty
.rc_original_make_medium_scenarios <- rc_make_medium_scenarios
.rc_original_run_microcompass <- rc_run_microcompass
.rc_original_run_regcompass <- rc_run_regcompass
.rc_original_run_regcompass_one_shot <- rc_run_regcompass_one_shot

.rc_clamp01 <- function(x) pmin(pmax(x, 0), 1)

.rc_absolute_activity_score <- function(X, half_saturation = 1) {
  X <- as.matrix(X)
  if (!is.numeric(half_saturation) || length(half_saturation) != 1L ||
      !is.finite(half_saturation) || half_saturation <= 0) {
    stop("`half_saturation` must be one positive finite number.", call. = FALSE)
  }
  observed <- is.finite(X)
  signal <- pmax(X, 0)
  score <- signal / (signal + half_saturation)
  score[observed & signal <= 0] <- 0
  score[!observed] <- NA_real_
  dimnames(score) <- dimnames(X)
  attr(score, "score_semantics") <- paste(
    "zero-preserving bounded support from non-negative normalized signal;",
    "not a probability or enzyme capacity"
  )
  score
}

.rc_relative_state_score <- function(X, min_scale = 0.05, z_clip = 6,
                                     tolerance = 1e-12) {
  X <- as.matrix(X)
  score <- rc_sigmoid(rc_gene_zscore(X, min_scale = min_scale,
                                     z_clip = z_clip))
  finite_range <- apply(X, 1L, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    diff(range(x))
  })
  noninformative <- !is.finite(finite_range) | finite_range <= tolerance
  if (any(noninformative)) score[noninformative, ] <- NA_real_
  dimnames(score) <- dimnames(X)
  attr(score, "score_semantics") <-
    "within_gene_relative_state_not_absolute_capacity"
  score
}

# The default is a zero-preserving absolute activity score. Relative state is
# retained as an explicit diagnostic and is NA for constant rows.
rc_gene_score <- function(X, min_scale = 0.05, z_clip = 6,
                          mode = c("absolute", "relative"),
                          half_saturation = getOption(
                            "RegCompassR.cpm_half_saturation", 1
                          )) {
  mode <- match.arg(mode)
  if (identical(mode, "absolute")) {
    return(.rc_absolute_activity_score(X, half_saturation))
  }
  .rc_relative_state_score(X, min_scale, z_clip)
}

.rc_weighted_gene_score <- function(X, weights, min_scale = 0.05,
                                    z_clip = 6,
                                    mode = c("absolute", "relative"),
                                    half_saturation = getOption(
                                      "RegCompassR.cpm_half_saturation", 1
                                    )) {
  X <- as.matrix(X)
  if (length(weights) != ncol(X) || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("`weights` must contain one positive finite value per column.",
         call. = FALSE)
  }
  mode <- match.arg(mode)
  if (identical(mode, "absolute")) {
    return(.rc_absolute_activity_score(X, half_saturation))
  }
  centers <- apply(X, 1L, .rc_weighted_quantile,
                   weights = weights, probs = 0.5)
  scales <- vapply(seq_len(nrow(X)), function(i) {
    mad_sigma <- .rc_weighted_quantile(
      abs(X[i, ] - centers[[i]]), weights, probs = 0.5
    ) * 1.4826
    quartiles <- .rc_weighted_quantile(
      X[i, ], weights, probs = c(0.25, 0.75)
    )
    max(mad_sigma, diff(quartiles) / 1.349, min_scale, na.rm = TRUE)
  }, numeric(1))
  z <- sweep(X, 1L, centers, "-")
  z <- sweep(z, 1L, scales, "/")
  z <- pmax(pmin(z, z_clip), -z_clip)
  score <- rc_sigmoid(z)
  finite_range <- apply(X, 1L, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    diff(range(x))
  })
  noninformative <- !is.finite(finite_range) | finite_range <= 1e-12
  if (any(noninformative)) score[noninformative, ] <- NA_real_
  dimnames(score) <- dimnames(X)
  attr(score, "score_semantics") <-
    "sample_balanced_within_gene_relative_state"
  score
}

.rc_capacity_diagnostics <- function(C_raw, q_values, relative,
                                     sample_balanced) {
  C_raw <- as.matrix(C_raw)
  n_finite <- rowSums(is.finite(C_raw))
  minimum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) min(x) else NA_real_
  })
  maximum <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) max(x) else NA_real_
  })
  data.frame(
    reaction_id = rownames(C_raw),
    stratum = if (sample_balanced) "global_sample_balanced" else "global",
    n = as.integer(n_finite),
    n_global = as.integer(n_finite),
    q_stratum = as.numeric(q_values),
    q_stratum_used = as.numeric(q_values),
    q_global = as.numeric(q_values),
    rho_n = 1,
    q_value = as.numeric(q_values),
    quantile_used = 0.95,
    n_finite = as.integer(n_finite),
    n_finite_global = as.integer(n_finite),
    low_n_flag = n_finite < 20L,
    all_missing_reaction_flag = n_finite == 0L,
    all_zero_reaction_flag = is.finite(maximum) & maximum <= 1e-12,
    constant_reaction_flag = is.finite(minimum) & is.finite(maximum) &
      abs(maximum - minimum) <= 1e-12,
    raw_out_of_unit_interval_flag = apply(C_raw, 1L, function(x) {
      any(is.finite(x) & (x < -1e-12 | x > 1 + 1e-12))
    }),
    relative_capacity_informative = rowSums(is.finite(relative)) > 0L,
    sample_balanced = sample_balanced,
    calibration_role = "diagnostic_only_not_lp_capacity",
    stringsAsFactors = FALSE
  )
}

# Q95 is retained as a within-reaction diagnostic only. For compatibility,
# C_rel now contains bounded absolute evidence used by the LP.
.rc_weighted_q95_calibrate <- function(C_raw, weights, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (length(weights) != ncol(C_raw) || any(!is.finite(weights)) ||
      any(weights <= 0)) {
    stop("`weights` must contain one positive finite value per capacity column.",
         call. = FALSE)
  }
  q_values <- apply(C_raw, 1L, .rc_weighted_quantile,
                    weights = weights, probs = 0.95)
  relative <- sweep(C_raw, 1L, q_values + eps, "/")
  relative <- .rc_clamp01(relative)
  noninformative <- !is.finite(q_values) | q_values <= eps |
    apply(C_raw, 1L, function(x) {
      x <- x[is.finite(x)]
      length(x) < 2L || diff(range(x)) <= eps
    })
  if (any(noninformative)) relative[noninformative, ] <- NA_real_

  absolute <- .rc_clamp01(C_raw)
  all_missing <- rowSums(is.finite(C_raw)) == 0L
  if (any(all_missing)) absolute[all_missing, ] <- NA_real_
  all_zero <- apply(C_raw, 1L, function(x) {
    x <- x[is.finite(x)]
    length(x) > 0L && max(x) <= eps
  })
  if (any(all_zero)) absolute[all_zero, ] <- 0

  list(
    C_rel = absolute,
    C_abs = absolute,
    C_within_reaction_relative = relative,
    Q = .rc_capacity_diagnostics(
      C_raw, q_values, relative, sample_balanced = TRUE
    )
  )
}

rc_q95_calibrate <- function(C_raw, eps = 1e-6, bootstrap = TRUE,
                             B = 500, BPPARAM = NULL, n0 = 80,
                             unit_meta = NULL, stratum_col = NULL) {
  C_raw <- as.matrix(C_raw)
  weights <- rep(1 / max(1L, ncol(C_raw)), ncol(C_raw))
  answer <- .rc_weighted_q95_calibrate(C_raw, weights, eps)
  answer$Q$sample_balanced <- FALSE
  answer$Q$stratum <- "global"
  if (isTRUE(bootstrap)) {
    answer$Q$q95_bootstrap <- NA_real_
    answer$Q$q95_ci_low <- NA_real_
    answer$Q$q95_ci_high <- NA_real_
    answer$Q$q95_ci_width <- NA_real_
    answer$Q$q95_unstable_flag <- NA
  }
  answer
}

.rc_gpr_tokenize <- function(gpr) {
  text <- gsub("([()])", " \\1 ", tolower(trimws(as.character(gpr))))
  tokens <- strsplit(gsub("\\s+", " ", text), " ", fixed = TRUE)[[1L]]
  tokens[nzchar(tokens)]
}

.rc_gpr_ast_to_dnf <- function(node, max_terms = 10000L) {
  if (identical(node$type, "gene")) return(list(node$value))
  child_terms <- lapply(node$children, .rc_gpr_ast_to_dnf,
                        max_terms = max_terms)
  if (identical(node$type, "or")) {
    answer <- unlist(child_terms, recursive = FALSE)
  } else {
    answer <- list(character())
    for (terms in child_terms) {
      answer <- unlist(lapply(answer, function(left) {
        lapply(terms, function(right) unique(c(left, right)))
      }), recursive = FALSE)
      if (length(answer) > max_terms) {
        stop("GPR expansion exceeded `max_terms`; use a structured long-table GPR.",
             call. = FALSE)
      }
    }
  }
  keys <- vapply(answer, function(x) {
    paste(sort(unique(x)), collapse = "\001")
  }, character(1))
  answer[!duplicated(keys)]
}

# Recursive Boolean parser with AND precedence over OR. Nested expressions are
# converted to the DNF list representation used by existing GPR functions.
rc_parse_gpr_simple <- function(gpr, max_terms = 10000L) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) return(list())
  tokens <- .rc_gpr_tokenize(gpr)
  position <- 1L
  current <- function() {
    if (position <= length(tokens)) tokens[[position]] else NA_character_
  }
  consume <- function(expected = NULL) {
    token <- current()
    if (!is.null(expected) && !identical(token, expected)) {
      found <- if (is.na(token)) "<end>" else token
      stop("Malformed GPR: expected `", expected, "` but found `",
           found, "`.", call. = FALSE)
    }
    position <<- position + 1L
    token
  }
  parse_primary <- NULL
  parse_and <- NULL
  parse_or <- NULL
  parse_primary <- function() {
    token <- current()
    if (is.na(token)) stop("Malformed GPR: unexpected end of rule.", call. = FALSE)
    if (identical(token, "(")) {
      consume("(")
      node <- parse_or()
      consume(")")
      return(node)
    }
    if (token %in% c("and", "or", ")")) {
      stop("Malformed GPR near token `", token, "`.", call. = FALSE)
    }
    consume()
    list(type = "gene", value = token)
  }
  parse_and <- function() {
    children <- list(parse_primary())
    while (identical(current(), "and")) {
      consume("and")
      children[[length(children) + 1L]] <- parse_primary()
    }
    if (length(children) == 1L) children[[1L]] else
      list(type = "and", children = children)
  }
  parse_or <- function() {
    children <- list(parse_and())
    while (identical(current(), "or")) {
      consume("or")
      children[[length(children) + 1L]] <- parse_and()
    }
    if (length(children) == 1L) children[[1L]] else
      list(type = "or", children = children)
  }
  ast <- parse_or()
  if (position <= length(tokens)) {
    stop("Malformed GPR: unexpected trailing token `", current(), "`.",
         call. = FALSE)
  }
  lapply(.rc_gpr_ast_to_dnf(ast, max_terms), function(x) {
    unique(x[nzchar(x)])
  })
}

# Human-GEM import now fails with reaction context instead of silently dropping
# parser failures.
rc_prepare_humangem_gpr_table <- function(repo_dir,
                                           gene_format = c("symbol", "ensembl")) {
  gene_format <- match.arg(gene_format)
  model_dir <- file.path(repo_dir, "model")
  genes_tsv <- file.path(model_dir, "genes.tsv")
  reactions_tsv <- file.path(model_dir, "reactions.tsv")
  model_yml <- file.path(model_dir, "Human-GEM.yml")
  missing <- c(genes_tsv, reactions_tsv, model_yml)[
    !file.exists(c(genes_tsv, reactions_tsv, model_yml))
  ]
  if (length(missing)) {
    stop("Missing Human-GEM model files: ",
         paste(basename(missing), collapse = ", "), call. = FALSE)
  }
  genes <- utils::read.delim(genes_tsv, stringsAsFactors = FALSE,
                             check.names = FALSE)
  if (!all(c("genes", "geneSymbols") %in% colnames(genes))) {
    stop("`genes.tsv` must contain `genes` and `geneSymbols` columns.",
         call. = FALSE)
  }
  gene_map <- stats::setNames(as.character(genes$geneSymbols),
                              as.character(genes$genes))
  bad <- is.na(gene_map) | !nzchar(gene_map)
  gene_map[bad] <- names(gene_map)[bad]
  reactions <- utils::read.delim(reactions_tsv, stringsAsFactors = FALSE,
                                 check.names = FALSE)
  if (!"rxns" %in% colnames(reactions)) {
    stop("`reactions.tsv` must contain an `rxns` column.", call. = FALSE)
  }
  reaction_rules <- rc_read_humangem_yml_rules(model_yml)
  reaction_rules <- reaction_rules[
    reaction_rules$reaction_id %in% reactions$rxns &
      nzchar(reaction_rules$gpr), , drop = FALSE
  ]
  rows <- lapply(seq_len(nrow(reaction_rules)), function(i) {
    reaction_id <- reaction_rules$reaction_id[[i]]
    rule <- reaction_rules$gpr[[i]]
    if (identical(gene_format, "symbol")) {
      rule <- rc_replace_humangem_gene_ids(rule, gene_map)
    }
    parsed <- tryCatch(
      rc_parse_gpr_simple(rule),
      error = function(error) {
        stop("Failed to parse Human-GEM GPR for reaction `", reaction_id,
             "`: ", conditionMessage(error), " Rule: ", rule,
             call. = FALSE)
      }
    )
    if (!length(parsed)) return(NULL)
    do.call(rbind, lapply(seq_along(parsed), function(group_index) {
      data.frame(
        reaction_id = reaction_id,
        and_group_id = group_index,
        gene = parsed[[group_index]],
        stringsAsFactors = FALSE
      )
    }))
  })
  nonempty <- rows[!vapply(rows, is.null, logical(1))]
  gpr_table <- if (length(nonempty)) do.call(rbind, nonempty) else
    data.frame(reaction_id = character(), and_group_id = integer(),
               gene = character())
  gpr_table$gene <- toupper(gpr_table$gene)
  gpr_table <- unique(gpr_table)
  rownames(gpr_table) <- NULL
  list(
    gpr_table = gpr_table,
    metabolic_genes = sort(unique(gpr_table$gene)),
    reaction_rules = reaction_rules,
    genes = genes,
    reactions = reactions,
    parser_diagnostics = data.frame(
      n_rules = nrow(reaction_rules),
      n_parsed = length(nonempty),
      n_failed = 0L,
      parser = "recursive_boolean_ast_to_dnf",
      stringsAsFactors = FALSE
    )
  )
}

# Canonical runs use bottleneck-preserving and annotation-count-neutral GPR
# defaults. Legacy heuristics remain available to direct callers.
rc_reaction_capacity <- function(
    gpr_list, gene_score,
    promiscuity_mode = c("sqrt", "linear", "none"),
    tau = 0.20,
    and_method = c("boltzmann", "min", "mean"),
    or_method = c("sum_sqrtK", "max", "prob_or", "sum"),
    BPPARAM = NULL) {
  if (isTRUE(getOption("RegCompassR.strict_gpr_defaults", FALSE))) {
    promiscuity_mode <- "none"
    and_method <- "min"
    or_method <- "max"
  }
  .rc_original_reaction_capacity(
    gpr_list, gene_score,
    promiscuity_mode = promiscuity_mode,
    tau = tau,
    and_method = and_method,
    or_method = or_method,
    BPPARAM = BPPARAM
  )
}

rc_map_meta_module_core_reactions <- local({
  original <- rc_map_meta_module_core_reactions
  function(gene_nodes, gpr_table) {
    if (!"and_group_id" %in% colnames(gpr_table)) {
      stop(
        "`gpr_table` must contain `and_group_id`; one group per gene would misclassify required subunits as isoenzymes.",
        call. = FALSE
      )
    }
    original(gene_nodes, gpr_table)
  }
})

.rc_case_insensitive_lookup <- function(ids) {
  keys <- toupper(trimws(as.character(ids)))
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  stats::setNames(as.character(ids)[keep], keys[keep])
}

# Signed TF x peak activity. A value of 0.5 is neutral, values below 0.5
# indicate active repression and values above 0.5 indicate active support.
.rc_pando_gene_confidence <- function(significant_edges, object, atac_assay,
                                      target_genes = NULL,
                                      rna_assay = "RNA") {
  units <- colnames(object)
  genes <- unique(toupper(as.character(target_genes)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes) && is.data.frame(significant_edges) &&
      "target" %in% colnames(significant_edges)) {
    genes <- unique(toupper(as.character(significant_edges$target)))
    genes <- genes[!is.na(genes) & nzchar(genes)]
  }
  confidence <- matrix(
    NA_real_, nrow = length(genes), ncol = length(units),
    dimnames = list(tolower(genes), units)
  )
  diagnostics <- data.frame(
    gene = genes,
    n_pando_edges = integer(length(genes)),
    n_positive_edges = integer(length(genes)),
    n_negative_edges = integer(length(genes)),
    n_unique_regions = integer(length(genes)),
    n_matched_regions = integer(length(genes)),
    n_unique_tfs = integer(length(genes)),
    n_matched_tfs = integer(length(genes)),
    matched_region_fraction = NA_real_,
    matched_tf_fraction = NA_real_,
    pando_supported = FALSE,
    confidence_source = "pando_signed_tf_peak_gene_regulatory_support",
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(significant_edges) || !nrow(significant_edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  required <- c("target", "region", "tf", "estimate")
  missing <- setdiff(required, colnames(significant_edges))
  if (length(missing)) {
    stop("Pando coefficient table lacks columns required for signed support: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  edges <- significant_edges
  edges$target <- toupper(trimws(as.character(edges$target)))
  edges$tf <- toupper(trimws(as.character(edges$tf)))
  edges$region <- trimws(as.character(edges$region))
  edges$estimate <- suppressWarnings(as.numeric(edges$estimate))
  edges <- edges[
    !is.na(edges$target) & nzchar(edges$target) &
      !is.na(edges$tf) & nzchar(edges$tf) &
      !is.na(edges$region) & nzchar(edges$region) &
      is.finite(edges$estimate) & edges$estimate != 0 &
      edges$target %in% genes, , drop = FALSE
  ]
  if (!nrow(edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  edges$.sign <- sign(edges$estimate)
  edges$.weight <- abs(edges$estimate)
  if ("rsq" %in% colnames(edges)) {
    rsq <- suppressWarnings(as.numeric(edges$rsq))
    quality <- sqrt(pmax(rsq, 0))
    quality[!is.finite(quality)] <- 0
    edges$.weight <- edges$.weight * quality
  }
  edges <- edges[is.finite(edges$.weight) & edges$.weight > 0, , drop = FALSE]
  if (!nrow(edges)) {
    return(list(gene_confidence = confidence, diagnostics = diagnostics))
  }
  atac <- .rc_pando_assay_data(object, atac_assay)
  rna <- .rc_pando_assay_data(object, rna_assay)
  peak_keys <- toupper(.rc_pando_region_key(rownames(atac)))
  keep_peak <- !is.na(peak_keys) & nzchar(peak_keys) & !duplicated(peak_keys)
  peak_lookup <- stats::setNames(rownames(atac)[keep_peak], peak_keys[keep_peak])
  edges$.peak_id <- unname(
    peak_lookup[toupper(.rc_pando_region_key(edges$region))]
  )
  tf_lookup <- .rc_case_insensitive_lookup(rownames(rna))
  edges$.tf_id <- unname(tf_lookup[edges$tf])

  for (gene in genes) {
    selected <- edges[edges$target == gene, , drop = FALSE]
    row <- diagnostics$gene == gene
    diagnostics$n_pando_edges[row] <- nrow(selected)
    diagnostics$n_positive_edges[row] <- sum(selected$.sign > 0)
    diagnostics$n_negative_edges[row] <- sum(selected$.sign < 0)
    diagnostics$n_unique_regions[row] <- length(unique(selected$region))
    diagnostics$n_unique_tfs[row] <- length(unique(selected$tf))
    matched <- selected[
      !is.na(selected$.peak_id) & nzchar(selected$.peak_id) &
        !is.na(selected$.tf_id) & nzchar(selected$.tf_id), , drop = FALSE
    ]
    diagnostics$n_matched_regions[row] <- length(unique(matched$.peak_id))
    diagnostics$n_matched_tfs[row] <- length(unique(matched$.tf_id))
    diagnostics$matched_region_fraction[row] <- if (
      diagnostics$n_unique_regions[row] > 0L
    ) diagnostics$n_matched_regions[row] / diagnostics$n_unique_regions[row] else
      NA_real_
    diagnostics$matched_tf_fraction[row] <- if (
      diagnostics$n_unique_tfs[row] > 0L
    ) diagnostics$n_matched_tfs[row] / diagnostics$n_unique_tfs[row] else
      NA_real_
    diagnostics$pando_supported[row] <- nrow(matched) > 0L
    if (!nrow(matched)) next
    peak_score <- rc_gene_score(
      as.matrix(atac[matched$.peak_id, units, drop = FALSE]),
      mode = "absolute",
      half_saturation = getOption("RegCompassR.atac_half_saturation", 1)
    )
    tf_score <- rc_gene_score(
      as.matrix(rna[matched$.tf_id, units, drop = FALSE]),
      mode = "absolute",
      half_saturation = getOption("RegCompassR.tf_half_saturation", 1)
    )
    edge_activity <- peak_score * tf_score
    weights <- matched$.weight / sum(matched$.weight)
    signed_activity <- as.numeric(crossprod(
      weights * matched$.sign, edge_activity
    ))
    confidence[tolower(gene), ] <- .rc_clamp01(
      0.5 + 0.5 * signed_activity
    )
  }
  list(gene_confidence = confidence, diagnostics = diagnostics)
}

.rc_pando_reaction_confidence <- function(meta_modules, pando_object, gem,
                                           atac_assay = "ATAC",
                                           rna_assay = "RNA") {
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  target_genes <- meta_modules$target_metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  gene <- .rc_pando_gene_confidence(
    meta_modules$tf_peak_gene_significant,
    pando_object,
    atac_assay = atac_assay,
    rna_assay = rna_assay,
    target_genes = target_genes
  )
  supported <- rowSums(is.finite(gene$gene_confidence)) > 0L
  gene_for_reaction <- gene$gene_confidence[supported, , drop = FALSE]
  reaction <- rc_reaction_confidence(
    parsed,
    gene_confidence = gene_for_reaction,
    unit_ids = colnames(pando_object)
  )
  list(
    gene_confidence = gene$gene_confidence,
    gene_confidence_diagnostics = gene$diagnostics,
    reaction_confidence = reaction,
    reaction_confidence_matrix = .rc_pando_reaction_confidence_matrix(
      reaction, names(parsed), colnames(pando_object)
    ),
    confidence_source = "pando_signed_tf_peak_gene_regulatory_support"
  )
}

.rc_signed_relation <- function(values) {
  values <- values[is.finite(values) & values != 0]
  if (!length(values)) return(NA_character_)
  if (all(values > 0)) return("concordant")
  if (all(values < 0)) return("discordant")
  "mixed"
}

# The original graph remains usable for components, but regulator, direction and
# sign metadata are retained for biological interpretation.
rc_project_metabolic_grn <- function(tf_peak_gene, metabolic_genes,
                                     top_k = 5L, min_shared_tfs = 1L,
                                     min_tf_jaccard = 0,
                                     max_targets_per_tf = 200L,
                                     include_direct_metabolic_tf = TRUE) {
  answer <- .rc_original_project_metabolic_grn(
    tf_peak_gene, metabolic_genes,
    top_k = top_k,
    min_shared_tfs = min_shared_tfs,
    min_tf_jaccard = min_tf_jaccard,
    max_targets_per_tf = max_targets_per_tf,
    include_direct_metabolic_tf = include_direct_metabolic_tf
  )
  edges <- answer$edges
  edges$regulator_set <- NA_character_
  edges$direct_regulator <- NA_character_
  edges$direct_target <- NA_character_
  edges$regulatory_relation <- NA_character_
  edges$signed_projection_weight <- NA_real_
  edges$direction_and_sign_preserved <- FALSE
  if (!nrow(edges) || !is.data.frame(tf_peak_gene) ||
      !all(c("sample_id", "tf", "target") %in% colnames(tf_peak_gene))) {
    answer$edges <- edges
    return(answer)
  }
  x <- tf_peak_gene
  x$sample_id <- as.character(x$sample_id)
  x$tf <- toupper(trimws(as.character(x$tf)))
  x$target <- toupper(trimws(as.character(x$target)))
  x$.estimate <- if ("estimate" %in% colnames(x))
    suppressWarnings(as.numeric(x$estimate)) else rep(NA_real_, nrow(x))
  x$.strength <- abs(x$.estimate)

  for (i in seq_len(nrow(edges))) {
    sample <- as.character(edges$sample_id[[i]])
    a <- toupper(as.character(edges$gene_a[[i]]))
    b <- toupper(as.character(edges$gene_b[[i]]))
    xs <- x[x$sample_id == sample, , drop = FALSE]
    direct <- xs[(xs$tf == a & xs$target == b) |
                   (xs$tf == b & xs$target == a), , drop = FALSE]
    if (nrow(direct)) {
      edges$direct_regulator[[i]] <- paste(unique(direct$tf), collapse = ";")
      edges$direct_target[[i]] <- paste(unique(direct$target), collapse = ";")
      edges$regulator_set[[i]] <- paste(unique(direct$tf), collapse = ";")
      edges$regulatory_relation[[i]] <- .rc_signed_relation(direct$.estimate)
      edges$signed_projection_weight[[i]] <- sum(direct$.estimate, na.rm = TRUE)
      edges$direction_and_sign_preserved[[i]] <- TRUE
      next
    }
    xa <- xs[xs$target == a, c("tf", ".estimate", ".strength"), drop = FALSE]
    xb <- xs[xs$target == b, c("tf", ".estimate", ".strength"), drop = FALSE]
    shared <- intersect(unique(xa$tf), unique(xb$tf))
    if (!length(shared)) next
    contributions <- vapply(shared, function(tf) {
      ea <- sum(xa$.estimate[xa$tf == tf], na.rm = TRUE)
      eb <- sum(xb$.estimate[xb$tf == tf], na.rm = TRUE)
      sa <- sum(xa$.strength[xa$tf == tf], na.rm = TRUE)
      sb <- sum(xb$.strength[xb$tf == tf], na.rm = TRUE)
      sign(ea) * sign(eb) * min(sa, sb)
    }, numeric(1))
    edges$regulator_set[[i]] <- paste(shared, collapse = ";")
    edges$regulatory_relation[[i]] <- .rc_signed_relation(contributions)
    edges$signed_projection_weight[[i]] <- sum(contributions, na.rm = TRUE)
    edges$direction_and_sign_preserved[[i]] <- TRUE
  }
  answer$edges <- edges
  answer
}

# Regulatory evidence is centered at 0.5. Missing/neutral regulation is neutral;
# only active repression contributes a non-negative confidence penalty.
.rc_compute_multiome_penalty_core <- function(
    C_rel, reaction_confidence, gpr_diagnostics = NULL,
    reaction_roles = NULL,
    weights = c(expr = 1.0, confidence = 0.5, missing = 1.0,
                gpr_missing = 0),
    eps = 1e-6, penalty_cap = 20,
    support_penalty = c(
      exchange = 0.05, demand = 20, sink = 20,
      artificial_support = 20, cofactor_recycle = 0.50, transport = 1.00
    ),
    missing_penalty = 5) {
  C <- as.matrix(C_rel)
  F_original <- rc_layer2_confidence_matrix(reaction_confidence, C)
  finite <- is.finite(F_original)
  F_original[finite] <- .rc_clamp01(F_original[finite])
  F_effective <- matrix(
    1, nrow = nrow(F_original), ncol = ncol(F_original),
    dimnames = dimnames(F_original)
  )
  F_effective[finite] <- pmin(1, pmax(2 * F_original[finite], eps))
  answer <- .rc_original_compute_multiome_penalty_core(
    C, F_effective,
    gpr_diagnostics = gpr_diagnostics,
    reaction_roles = reaction_roles,
    weights = weights,
    eps = eps,
    penalty_cap = penalty_cap,
    support_penalty = support_penalty,
    missing_penalty = missing_penalty
  )
  answer$components$reaction_regulatory_support <- F_original
  answer$components$reaction_confidence_effective <- F_effective
  answer$components$missing_regulatory_support_flag <- !finite
  answer$evidence_policy <- paste(
    "zero-preserving absolute RNA evidence with signed Pando modulation;",
    "0.5 is regulatory neutral, missing regulation is neutral, and only",
    "active repression contributes a confidence penalty"
  )
  answer$penalty_formula <- paste(
    "w_expr*-log(C_abs) +",
    "w_conf*-log(min(1,max(2*regulatory_support,eps))) +",
    "w_missing*missing_expression + w_gpr*gpr_missing"
  )
  answer
}

# Structural support cannot silently replace biological evidence.
rc_layer2_penalty <- function(C_rel, Conf, epsilon = 1e-6,
                              epsilon_C = 1e-3, epsilon_Conf = 1e-3,
                              penalty_cap = 20,
                              support_reactions = character(),
                              support_penalty = 0.05,
                              allow_structural_support_override = FALSE) {
  if (length(support_reactions) && !isTRUE(allow_structural_support_override)) {
    stop(
      "FASTCORE/support membership is structural, not biological evidence. Set `allow_structural_support_override = TRUE` only for sensitivity analysis.",
      call. = FALSE
    )
  }
  .rc_original_layer2_penalty(
    C_rel, Conf,
    epsilon = epsilon,
    epsilon_C = epsilon_C,
    epsilon_Conf = epsilon_Conf,
    penalty_cap = penalty_cap,
    support_reactions = support_reactions,
    support_penalty = support_penalty
  )
}

# Stable within-target display rank. It is not a probability. Constant targets
# are non-informative and remain NA.
rc_compass_score_from_penalty <- function(P, feasible, epsilon = 1e-6,
                                          method = c("ecdf", "mad_sigmoid"),
                                          variation_tolerance = 1e-8) {
  method <- match.arg(method)
  P <- as.matrix(P)
  feasible <- as.matrix(feasible)
  score <- matrix(NA_real_, nrow(P), ncol(P), dimnames = dimnames(P))
  noninformative <- logical(nrow(P))
  for (i in seq_len(nrow(P))) {
    index <- feasible[i, ] & is.finite(P[i, ])
    x <- P[i, index]
    if (length(x) < 2L || diff(range(x)) <= variation_tolerance) {
      noninformative[[i]] <- TRUE
      next
    }
    if (identical(method, "ecdf")) {
      ranks <- rank(x, ties.method = "average")
      score[i, index] <- 1 - (ranks - 1) / (length(x) - 1)
    } else {
      center <- stats::median(x)
      scale <- max(stats::mad(x, constant = 1.4826),
                   stats::IQR(x) / 1.349, variation_tolerance)
      score[i, index] <- rc_sigmoid((center - x) / scale)
    }
  }
  attr(score, "score_semantics") <- if (identical(method, "ecdf")) {
    "within_target_relative_penalty_rank_not_probability"
  } else {
    "within_target_robust_penalty_transform_not_probability"
  }
  attr(score, "noninformative_target") <- stats::setNames(
    noninformative, rownames(P)
  )
  score
}

# Medium backgrounds use only the current explicit scenario identifiers.
rc_make_medium_scenarios <- function(
    gem,
    scenario = "permissive_all_exchange",
    custom_medium = NULL,
    uptake_scale = c(
      permissive_all_exchange = 1,
      normal_human_plasma = 1, rpmi1640 = 1, minimal = 0.1,
      low_glucose = 0.1,
      low_glutamine = 0.1, high_lactate = 1
    ),
    condition_col = NULL,
    exchange_roles = c("exchange"),
    condition = condition_col) {
  choices <- c(
    "permissive_all_exchange", "normal_human_plasma", "rpmi1640", "minimal",
    "low_glucose", "low_glutamine",
    "high_lactate", "custom"
  )
  scenario <- match.arg(scenario, choices = choices, several.ok = TRUE)
  if ("custom" %in% scenario) {
    if (is.null(custom_medium)) {
      stop("`custom_medium` is required when `scenario` includes 'custom'.",
           call. = FALSE)
    }
    required <- c("medium_scenario_id", "exchange_reaction_id",
                  "lb", "ub", "available")
    missing <- setdiff(required, colnames(custom_medium))
    if (length(missing)) {
      stop("`custom_medium` missing columns: ", paste(missing, collapse = ", "),
           call. = FALSE)
    }
  }
  baseline <- .rc_original_make_medium_scenarios(
    gem, scenario = "normal_human_plasma", uptake_scale = 1,
    condition_col = condition_col,
    exchange_roles = exchange_roles,
    condition = condition
  )
  annotation <- tolower(paste(
    baseline$exchange_reaction_id,
    if ("metabolite_id" %in% colnames(baseline)) baseline$metabolite_id else ""
  ))
  names(annotation) <- baseline$exchange_reaction_id
  scale_value <- function(name, default) {
    if (!is.null(names(uptake_scale)) && name %in% names(uptake_scale)) {
      as.numeric(uptake_scale[[name]])
    } else if (length(uptake_scale) == 1L) {
      as.numeric(uptake_scale[[1L]])
    } else default
  }
  build_one <- function(name) {
    if (identical(name, "custom")) return(NULL)
    out <- baseline
    out$medium_scenario_id <- name
    base_scale <- if (identical(name, "minimal"))
      scale_value("minimal", 0.1) else
      scale_value("permissive_all_exchange", 1)
    reaction_scale <- rep(base_scale, nrow(out))
    target_pattern <- switch(
      name,
      low_glucose = "glucose|d[- ]?glucose|glc",
      low_glutamine = "glutamine|gln",
      high_lactate = "lactate|lac",
      NULL
    )
    target <- rep(FALSE, nrow(out))
    if (!is.null(target_pattern)) {
      target <- grepl(target_pattern,
                      annotation[out$exchange_reaction_id], perl = TRUE)
      if (!any(target)) {
        stop("Scenario `", name,
             "` could not identify its target exchange; supply `custom_medium`.",
             call. = FALSE)
      }
      reaction_scale[target] <- scale_value(name, 1)
    }
    out$lb <- -10 * reaction_scale
    out$ub <- 1000
    out$available <- TRUE
    out$target_exchange_flag <- target
    out$evidence_source <- if (identical(name, "permissive_all_exchange")) {
      "explicit_permissive_all_exchange_technical_baseline"
    } else if (name %in% c("normal_human_plasma", "rpmi1640")) {
      "current_named_medium_definition"
    } else {
      "permissive_baseline_sensitivity_scenario"
    }
    out$assumption_level <- if (identical(name, "permissive_all_exchange")) {
      "technical_upper_bound"
    } else if (name %in% c("normal_human_plasma", "rpmi1640")) {
      "named_medium_background"
    } else {
      "sensitivity_only"
    }
    out$concentration_used_for_rate_bound <- FALSE
    out$rate_bound_source <- "unmeasured_assumption"
    out
  }
  built <- lapply(setdiff(scenario, "custom"), build_one)
  built <- built[!vapply(built, is.null, logical(1))]
  output <- if (length(built)) do.call(rbind, built) else NULL
  if ("custom" %in% scenario) {
    custom <- custom_medium
    optional <- c(
      "metabolite_id", "condition", "evidence_source", "assumption_level",
      "target_exchange_flag", "concentration_used_for_rate_bound",
      "rate_bound_source"
    )
    for (name in setdiff(optional, colnames(custom))) custom[[name]] <- NA
    if (is.null(output)) {
      output <- custom
    } else {
      columns <- union(colnames(output), colnames(custom))
      for (name in setdiff(columns, colnames(output))) output[[name]] <- NA
      for (name in setdiff(columns, colnames(custom))) custom[[name]] <- NA
      output <- rbind(output[, columns, drop = FALSE],
                      custom[, columns, drop = FALSE])
    }
  }
  rownames(output) <- NULL
  output
}

# Sample x cell type is the inference default. Metacell scores remain available
# for exploratory within-sample heterogeneity.
rc_run_microcompass <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("sample_celltype", "metacell"),
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type", model_params = list(),
    penalty_weights = c(expr = 1.0, confidence = 0.5, missing = 1.0),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    BPPARAM = NULL) {
  forced_unit <- getOption("RegCompassR.inference_unit", NULL)
  if (!is.null(forced_unit)) unit <- forced_unit
  unit <- match.arg(unit)
  if (identical(unit, "metacell")) {
    warning(
      "Metacell-level scores are exploratory within-sample observations. Use `unit = 'sample_celltype'` for biological-replicate inference.",
      call. = FALSE
    )
  }
  answer <- .rc_original_run_microcompass(
    layer1, gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = mode,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    penalty_weights = penalty_weights,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM
  )
  answer$relative_penalty_rank <- answer$score
  answer$score_semantics <- attr(answer$score, "score_semantics") %||%
    "within_target_relative_penalty_rank_not_probability"
  answer$noninformative_target <- attr(answer$score,
                                        "noninformative_target")
  answer$primary_output <- "penalty"
  answer$primary_output_semantics <-
    "minimum evidence-discordance penalty; lower means stronger support"
  answer$params$inference_unit <- unit
  answer
}

rc_run_regcompass <- function(
    object, gem, outdir, pfm, genome,
    fragment_files = NULL,
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "cell_type",
    rna_assay = "RNA",
    atac_assay = "ATAC",
    model_mode = c("meta_module_gem", "full_gem"),
    medium_scenarios = NULL,
    metacell_args = list(),
    layer1_args = list(),
    pando_args = list(),
    layer2_args = list(),
    upstream_workers = NULL,
    layer2_workers = NULL,
    parallel_backend = c("auto", "serial", "snow", "multicore"),
    strict_biological_defaults = TRUE,
    inference_unit = c("sample_celltype", "metacell")) {
  inference_unit <- match.arg(inference_unit)
  if (isTRUE(strict_biological_defaults)) {
    if (is.null(layer1_args$promiscuity_mode)) {
      layer1_args$promiscuity_mode <- "none"
    }
    if (is.null(layer1_args$and_method)) layer1_args$and_method <- "min"
  }
  previous <- options(
    RegCompassR.strict_gpr_defaults = isTRUE(strict_biological_defaults),
    RegCompassR.inference_unit = inference_unit
  )
  on.exit(options(previous), add = TRUE)
  answer <- .rc_original_run_regcompass(
    object, gem, outdir, pfm, genome,
    fragment_files = fragment_files,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    model_mode = model_mode,
    medium_scenarios = medium_scenarios,
    metacell_args = metacell_args,
    layer1_args = layer1_args,
    pando_args = pando_args,
    layer2_args = layer2_args,
    upstream_workers = upstream_workers,
    layer2_workers = layer2_workers,
    parallel_backend = parallel_backend
  )
  answer$schema_version <- "regcompass_v4_absolute_sample_inference"
  answer$layer1$C_abs <- answer$layer1$C_rel
  calibration <- .rc_weighted_q95_calibrate(
    answer$layer1$C_raw,
    answer$layer1$sample_balance_weights
  )
  answer$layer1$C_within_reaction_relative <-
    calibration$C_within_reaction_relative
  answer$layer1$gene_activity_absolute <- answer$layer1$global_gene_score
  answer$layer1$gene_state_relative <- .rc_weighted_gene_score(
    answer$layer1$rna_metacell_logcpm,
    answer$layer1$sample_balance_weights,
    mode = "relative"
  )
  answer$layer1$capacity_calibration_scope <-
    "zero_preserving_absolute_activity_with_q95_diagnostic_only"
  if (isTRUE(strict_biological_defaults)) {
    answer$layer1$capacity_params$promiscuity_mode <- "none"
    answer$layer1$capacity_params$and_method <- "min"
    answer$layer1$capacity_params$or_method <- "max"
  }
  answer$layer1$reaction_confidence_source <-
    "pando_signed_tf_peak_gene_regulatory_support"
  answer$params$penalty_unit <- inference_unit
  answer$params$strict_biological_defaults <- strict_biological_defaults
  answer$params$gpr_and_method <- if (strict_biological_defaults) "min" else
    layer1_args$and_method %||% "boltzmann"
  answer$params$gpr_or_method <- if (strict_biological_defaults) "max" else
    "sum_sqrtK"
  answer$params$promiscuity_mode <- if (strict_biological_defaults) "none" else
    layer1_args$promiscuity_mode %||% "sqrt"
  answer$params$sample_balanced_q95 <- FALSE
  answer$params$q95_role <- "diagnostic_only"
  saveRDS(answer, file.path(outdir, "regcompass_global_metacell_result.rds"))
  saveRDS(answer, file.path(outdir, "regcompass_result.rds"))
  answer
}

rc_run_regcompass_one_shot <- function(
    object, outdir, pfm, genome,
    fragment_files = NULL,
    gem = NULL,
    humangem_version = "2.0.0",
    medium_scenario = "permissive_all_exchange",
    medium_scenarios = NULL,
    ...) {
  .rc_original_run_regcompass_one_shot(
    object, outdir, pfm, genome,
    fragment_files = fragment_files,
    gem = gem,
    humangem_version = humangem_version,
    medium_scenario = medium_scenario,
    medium_scenarios = medium_scenarios,
    ...
  )
}
