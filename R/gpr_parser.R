#' Parse a simple GPR rule
#'
#' Supports flat GPR forms such as `(g1 and g2) or g3`, `g1 or g2`,
#' `g1 and g2`, and `g1`. Complex nested mixtures such as
#' `g1 and (g2 or g3)` stop rather than being silently mis-parsed.
rc_parse_gpr_simple <- function(gpr) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) return(list())
  gpr <- gsub("\\s+", " ", trimws(tolower(as.character(gpr))))
  if (grepl("and\\s*\\(", gpr) || grepl("\\)\\s*and", gpr)) {
    stop(
      "Complex nested AND/OR GPR formulas are not supported; provide a long-table GPR.",
      call. = FALSE
    )
  }
  flat <- gsub("\\([^()]+\\)", "X", gpr)
  if (grepl("[()]", flat)) {
    stop(
      "Complex or nested GPR formulas are not supported; provide a long-table GPR.",
      call. = FALSE
    )
  }
  gpr <- gsub("[()]", "", gpr)
  or_parts <- strsplit(gpr, "\\s+or\\s+", perl = TRUE)[[1L]]
  if (!length(or_parts) || any(!nzchar(trimws(or_parts)))) {
    stop("Malformed GPR rule contains an empty OR group.", call. = FALSE)
  }
  groups <- lapply(or_parts, function(group) {
    genes <- trimws(strsplit(group, "\\s+and\\s+", perl = TRUE)[[1L]])
    if (!length(genes) || any(!nzchar(genes))) {
      stop("Malformed GPR rule contains an empty AND subunit.", call. = FALSE)
    }
    unique(genes)
  })
  groups
}

#' Parse a reaction GPR table
#'
#' Accepts either `reaction_id + gpr` or a pre-parsed long table with
#' `reaction_id + and_group_id + gene`.
rc_parse_gpr_table <- function(gpr_table) {
  if (!is.data.frame(gpr_table)) {
    stop("`gpr_table` must be a data.frame.", call. = FALSE)
  }
  if (all(c("reaction_id", "and_group_id", "gene") %in% colnames(gpr_table))) {
    reaction_id <- trimws(as.character(gpr_table$reaction_id))
    and_group_id <- trimws(as.character(gpr_table$and_group_id))
    gene <- trimws(tolower(as.character(gpr_table$gene)))
    bad <- is.na(reaction_id) | !nzchar(reaction_id) |
      is.na(and_group_id) | !nzchar(and_group_id) |
      is.na(gene) | !nzchar(gene)
    if (any(bad)) {
      stop(
        "Long-table GPR reaction_id, and_group_id, and gene values must be non-missing and non-empty.",
        call. = FALSE
      )
    }
    normalized <- unique(data.frame(
      reaction_id = reaction_id,
      and_group_id = and_group_id,
      gene = gene,
      stringsAsFactors = FALSE
    ))
    reaction_levels <- unique(normalized$reaction_id)
    split_rxn <- split(
      normalized,
      factor(normalized$reaction_id, levels = reaction_levels),
      drop = TRUE
    )
    return(lapply(split_rxn, function(df) {
      group_levels <- unique(df$and_group_id)
      split_gene <- split(
        df$gene,
        factor(df$and_group_id, levels = group_levels),
        drop = TRUE
      )
      unname(lapply(split_gene, unique))
    }))
  }

  missing_cols <- setdiff(c("reaction_id", "gpr"), colnames(gpr_table))
  if (length(missing_cols) > 0L) {
    stop(
      "`gpr_table` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  reaction_id <- trimws(as.character(gpr_table$reaction_id))
  if (anyNA(reaction_id) || any(!nzchar(reaction_id))) {
    stop("`reaction_id` values must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(reaction_id)) {
    duplicated_ids <- unique(reaction_id[duplicated(reaction_id)])
    stop(
      "Simple GPR tables must contain one row per reaction_id; duplicated IDs: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  parsed <- lapply(gpr_table$gpr, rc_parse_gpr_simple)
  names(parsed) <- reaction_id
  parsed
}

#' Compute gene promiscuity weights from parsed GPR rules
rc_promiscuity_weight <- function(gpr_list, mode = c("sqrt", "linear", "none")) {
  mode <- match.arg(mode)
  if (!is.list(gpr_list)) stop("`gpr_list` must be a list.", call. = FALSE)
  reaction_genes <- lapply(gpr_list, function(rule) {
    unique(tolower(trimws(as.character(unlist(rule, use.names = FALSE)))))
  })
  reaction_genes <- lapply(reaction_genes, function(x) x[!is.na(x) & nzchar(x)])
  genes <- unique(unlist(reaction_genes, use.names = FALSE))
  if (!length(genes)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  gene_rxn <- table(unlist(reaction_genes, use.names = FALSE))
  n_rxn <- as.numeric(gene_rxn)
  weights <- switch(
    mode,
    none = rep(1, length(n_rxn)),
    sqrt = 1 / sqrt(n_rxn),
    linear = 1 / n_rxn
  )
  stats::setNames(weights, names(gene_rxn))
}
