# Focused corrections for core RegCompass workflow semantics.

.rc_expand_meta_module_reactions_uncorrected <- rc_expand_meta_module_reactions
.rc_build_meta_module_gem_uncorrected <- rc_build_meta_module_gem
.rc_compute_multiome_penalty_uncorrected <- rc_compute_multiome_penalty

.rc_hard_core_rows <- function(core_reactions) {
  if (is.null(core_reactions) || !is.data.frame(core_reactions)) {
    return(core_reactions)
  }
  if ("is_core" %in% colnames(core_reactions)) {
    return(core_reactions[core_reactions$is_core %in% TRUE, , drop = FALSE])
  }
  core_reactions
}

rc_build_metacell_metadata <- function(membership,
                                       metacell_id_col = "metacell_id",
                                       cell_id_col = "cell_id") {
  if (!is.data.frame(membership)) {
    stop("`membership` must be a data.frame.", call. = FALSE)
  }
  missing <- setdiff(
    c(metacell_id_col, cell_id_col),
    colnames(membership)
  )
  if (length(missing)) {
    stop("`membership` is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  metacell_id <- trimws(as.character(membership[[metacell_id_col]]))
  keep <- !is.na(metacell_id) & nzchar(metacell_id)
  x <- membership[keep, , drop = FALSE]
  x[[metacell_id_col]] <- metacell_id[keep]
  if (!nrow(x)) {
    out <- x[, setdiff(colnames(x), cell_id_col), drop = FALSE]
    out$n_cells <- integer()
    return(out)
  }

  strict_columns <- intersect(
    c("sample_id", "condition", "cell_type"),
    colnames(x)
  )
  split_rows <- split(seq_len(nrow(x)), x[[metacell_id_col]])
  for (metacell in names(split_rows)) {
    rows <- split_rows[[metacell]]
    for (column in strict_columns) {
      values <- trimws(as.character(x[[column]][rows]))
      values <- unique(values[!is.na(values) & nzchar(values)])
      if (length(values) != 1L || anyNA(x[[column]][rows]) ||
          any(!nzchar(trimws(as.character(x[[column]][rows]))))) {
        stop(
          "Metacell `", metacell, "` mixes metadata or contains missing values in `",
          column, "`.",
          call. = FALSE
        )
      }
    }
  }

  columns <- setdiff(colnames(x), cell_id_col)
  out <- x[!duplicated(x[[metacell_id_col]]), columns, drop = FALSE]
  out$n_cells <- as.integer(vapply(
    as.character(out[[metacell_id_col]]),
    function(id) length(split_rows[[id]]),
    integer(1)
  ))
  rownames(out) <- NULL
  out
}

rc_q95_shrink <- function(C_raw, unit_meta = NULL, stratum_col = NULL,
                          q = 0.95, n0 = 80, eps = 1e-6) {
  C_raw <- as.matrix(C_raw)
  if (is.null(rownames(C_raw)) || is.null(colnames(C_raw))) {
    stop("`C_raw` must have reaction rownames and unit colnames.", call. = FALSE)
  }
  global_q <- apply(C_raw, 1, rc_safe_quantile, probs = q)
  global_n <- rowSums(is.finite(C_raw))

  if (!is.null(stratum_col)) {
    if (is.null(unit_meta) || !"pool_id" %in% colnames(unit_meta) ||
        !stratum_col %in% colnames(unit_meta)) {
      stop("`unit_meta` with `pool_id` and `stratum_col` is required.",
           call. = FALSE)
    }
    unit_meta <- unit_meta[match(colnames(C_raw), unit_meta$pool_id), , drop = FALSE]
    if (anyNA(unit_meta$pool_id)) {
      stop("`unit_meta` is missing metadata for some capacity columns.",
           call. = FALSE)
    }
    strata <- trimws(as.character(unit_meta[[stratum_col]]))
    if (anyNA(strata) || any(!nzchar(strata))) {
      stop("Q95 strata must be non-missing and non-empty.", call. = FALSE)
    }
    if (length(unique(strata)) < 2L) {
      warning(
        "Only one stratum detected; Q95 stratified calibration degenerates to global calibration.",
        call. = FALSE
      )
    }
  } else {
    strata <- rep("global", ncol(C_raw))
  }

  diagnostics <- list()
  C_rel <- C_raw
  index <- 1L
  for (stratum in unique(strata)) {
    pools <- which(strata == stratum)
    n <- rowSums(is.finite(C_raw[, pools, drop = FALSE]))
    rho <- n / (n + n0)
    q_stratum <- apply(
      C_raw[, pools, drop = FALSE],
      1,
      rc_safe_quantile,
      probs = q
    )
    q_used <- ifelse(is.finite(q_stratum), q_stratum, global_q)
    q_shrink <- rho * q_used + (1 - rho) * global_q
    C_rel[, pools] <- sweep(
      C_raw[, pools, drop = FALSE],
      1,
      q_shrink + eps,
      "/"
    )
    diagnostics[[index]] <- data.frame(
      reaction_id = rownames(C_raw),
      stratum = stratum,
      n = as.integer(n),
      n_global = as.integer(global_n),
      q_stratum = as.numeric(q_stratum),
      q_stratum_used = as.numeric(q_used),
      q_global = as.numeric(global_q),
      rho_n = as.numeric(rho),
      q_shrink = as.numeric(q_shrink),
      q95_power_class = factor(
        ifelse(
          n < 5L,
          "very_low",
          ifelse(
            n < 20L,
            "low",
            ifelse(n < 100L, "moderate",
                   ifelse(n < 400L, "adequate", "high"))
          )
        ),
        levels = c("very_low", "low", "moderate", "adequate", "high"),
        ordered = TRUE
      ),
      stringsAsFactors = FALSE
    )
    index <- index + 1L
  }

  C_rel[C_rel > 1] <- 1
  all_missing <- global_n == 0L
  if (any(all_missing)) C_rel[all_missing, ] <- NA_real_
  Q <- do.call(rbind, diagnostics)
  Q$all_missing_reaction_flag <- Q$n_global == 0L
  Q$stratum_missing_reaction_flag <- Q$n == 0L
  list(C_rel = C_rel, Q = Q)
}

rc_layer2_has_gpr <- function(meta, C_rel = NULL, Conf = NULL) {
  if (!is.data.frame(meta) || !"reaction_id" %in% colnames(meta)) {
    stop("`meta` must contain `reaction_id`.", call. = FALSE)
  }
  gpr_columns <- intersect(
    c("gpr", "grRule", "gene_reaction_rule",
      "gene_reaction_rule_string", "genes"),
    colnames(meta)
  )
  has_metadata_gpr <- if (length(gpr_columns)) {
    apply(meta[, gpr_columns, drop = FALSE], 1, function(value) {
      any(!is.na(value) & nzchar(trimws(as.character(value))))
    })
  } else {
    rep(FALSE, nrow(meta))
  }

  has_evidence <- rep(FALSE, nrow(meta))
  if (!is.null(C_rel) && !is.null(Conf)) {
    C <- as.matrix(C_rel)
    F <- as.matrix(Conf)
    if (!is.null(rownames(C)) && !is.null(rownames(F))) {
      reaction_id <- as.character(meta$reaction_id)
      present <- reaction_id %in% rownames(C) & reaction_id %in% rownames(F)
      if (any(present)) {
        selected <- reaction_id[present]
        has_evidence[present] <- rowSums(
          is.finite(C[selected, , drop = FALSE]) |
            is.finite(F[selected, , drop = FALSE])
        ) > 0
      }
    }
  }
  has_metadata_gpr | has_evidence
}

rc_map_meta_module_core_reactions <- function(gene_nodes, gpr_table) {
  if (!is.data.frame(gene_nodes) ||
      !all(c("sample_id", "gene", "module_id") %in% colnames(gene_nodes))) {
    stop("`gene_nodes` must contain sample_id, gene and module_id.",
         call. = FALSE)
  }
  if (!is.data.frame(gpr_table) ||
      !all(c("reaction_id", "gene") %in% colnames(gpr_table))) {
    stop("`gpr_table` must contain reaction_id and gene.", call. = FALSE)
  }

  nodes <- unique(gene_nodes[, c("sample_id", "module_id", "gene"), drop = FALSE])
  nodes$sample_id <- as.character(nodes$sample_id)
  nodes$module_id <- as.character(nodes$module_id)
  nodes$gene <- toupper(trimws(as.character(nodes$gene)))
  nodes <- nodes[!is.na(nodes$gene) & nzchar(nodes$gene), , drop = FALSE]

  gpr <- gpr_table
  gpr$reaction_id <- trimws(as.character(gpr$reaction_id))
  gpr$gene <- toupper(trimws(as.character(gpr$gene)))
  gpr <- gpr[
    !is.na(gpr$reaction_id) & nzchar(gpr$reaction_id) &
      !is.na(gpr$gene) & nzchar(gpr$gene),
    , drop = FALSE
  ]
  if (!"and_group_id" %in% colnames(gpr)) {
    gpr$and_group_id <- ave(
      seq_len(nrow(gpr)),
      gpr$reaction_id,
      FUN = seq_along
    )
  }
  gpr$and_group_id <- as.character(gpr$and_group_id)
  gpr <- unique(gpr[, c("reaction_id", "and_group_id", "gene"), drop = FALSE])

  empty <- data.frame(
    sample_id = character(),
    module_id = character(),
    gene = character(),
    reaction_id = character(),
    and_group_id = character(),
    required_genes = character(),
    matched_genes = character(),
    missing_genes = character(),
    group_complete = logical(),
    is_core = logical(),
    is_partial_candidate = logical(),
    inclusion_stage = character(),
    stringsAsFactors = FALSE
  )
  if (!nrow(nodes) || !nrow(gpr)) return(empty)

  module_rows <- split(
    seq_len(nrow(nodes)),
    paste(nodes$sample_id, nodes$module_id, sep = "\001")
  )
  output <- list()
  output_index <- 0L

  for (module_key in names(module_rows)) {
    module_index <- module_rows[[module_key]]
    module_genes <- unique(nodes$gene[module_index])
    sample_id <- nodes$sample_id[module_index[[1L]]]
    module_id <- nodes$module_id[module_index[[1L]]]
    candidate_reactions <- unique(gpr$reaction_id[gpr$gene %in% module_genes])

    for (reaction_id in candidate_reactions) {
      reaction_gpr <- gpr[gpr$reaction_id == reaction_id, , drop = FALSE]
      group_rows <- split(seq_len(nrow(reaction_gpr)), reaction_gpr$and_group_id)
      group_complete <- vapply(group_rows, function(rows) {
        all(unique(reaction_gpr$gene[rows]) %in% module_genes)
      }, logical(1))
      reaction_core <- any(group_complete)

      for (group_id in names(group_rows)) {
        rows <- group_rows[[group_id]]
        required <- sort(unique(reaction_gpr$gene[rows]))
        matched <- intersect(required, module_genes)
        if (!length(matched)) next
        missing <- setdiff(required, module_genes)
        complete <- !length(missing)

        for (gene in matched) {
          output_index <- output_index + 1L
          output[[output_index]] <- data.frame(
            sample_id = sample_id,
            module_id = module_id,
            gene = gene,
            reaction_id = reaction_id,
            and_group_id = group_id,
            required_genes = paste(required, collapse = ";"),
            matched_genes = paste(sort(matched), collapse = ";"),
            missing_genes = paste(sort(missing), collapse = ";"),
            group_complete = complete,
            is_core = reaction_core,
            is_partial_candidate = !reaction_core,
            inclusion_stage = if (reaction_core) {
              "core_complete_gpr"
            } else {
              "partial_gpr_candidate"
            },
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(output)) return(empty)
  answer <- unique(do.call(rbind, output))
  rownames(answer) <- NULL
  answer
}

rc_expand_meta_module_reactions <- function(gem, core_reactions,
                                             subsystem_table = NULL,
                                             expansion_mode = c(
                                               "ordered_once", "fixed_point"
                                             ),
                                             max_iterations = 10L) {
  hard_core_reactions <- .rc_hard_core_rows(core_reactions)
  answer <- .rc_expand_meta_module_reactions_uncorrected(
    gem = gem,
    core_reactions = hard_core_reactions,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )
  if (!"is_core" %in% colnames(core_reactions)) return(answer)

  key <- function(data) paste(
    as.character(data$sample_id),
    as.character(data$module_id),
    as.character(data$reaction_id),
    sep = "\001"
  )
  candidate_keys <- unique(key(core_reactions))
  hard_rows <- .rc_hard_core_rows(core_reactions)
  hard_keys <- unique(key(hard_rows))
  membership_keys <- key(answer$reaction_membership)

  answer$reaction_membership$is_core <- membership_keys %in% hard_keys
  partial_anchor <- membership_keys %in%
    setdiff(candidate_keys, hard_keys)
  if (any(partial_anchor)) {
    answer$reaction_membership <- answer$reaction_membership[
      !partial_anchor,
      , drop = FALSE
    ]
    membership_keys <- key(answer$reaction_membership)
  }

  if (is.data.frame(answer$summary) && nrow(answer$summary)) {
    for (i in seq_len(nrow(answer$summary))) {
      selected <- as.character(core_reactions$sample_id) ==
        as.character(answer$summary$sample_id[[i]]) &
        as.character(core_reactions$module_id) ==
        as.character(answer$summary$module_id[[i]]) &
        core_reactions$is_core %in% TRUE
      answer$summary$n_core_reactions[[i]] <- length(unique(
        as.character(core_reactions$reaction_id[selected])
      ))
      gene_selected <- selected
      if ("group_complete" %in% colnames(core_reactions)) {
        gene_selected <- gene_selected & core_reactions$group_complete %in% TRUE
      }
      answer$summary$n_core_genes[[i]] <- length(unique(
        as.character(core_reactions$gene[gene_selected])
      ))
    }
  }
  answer
}

rc_build_meta_module_gem <- function(gem, reaction_membership,
                                     core_reactions = NULL, ...) {
  if (!is.null(core_reactions)) {
    core_reactions <- .rc_hard_core_rows(core_reactions)
  }
  .rc_build_meta_module_gem_uncorrected(
    gem = gem,
    reaction_membership = reaction_membership,
    core_reactions = core_reactions,
    ...
  )
}

rc_compute_multiome_penalty <- function(...) {
  answer <- .rc_compute_multiome_penalty_uncorrected(...)
  answer$evidence_policy <- "penalty_only"
  answer$evidence_description <- paste(
    "Multiome evidence modifies the LP objective penalty only;",
    "it does not directly change stoichiometry or reaction bounds."
  )
  answer
}
