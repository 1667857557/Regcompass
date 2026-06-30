#' Parse a simple GPR rule
#'
#' Supports v0.3 GPR forms such as `(g1 and g2) or g3`, `g1 or g2`,
#' `g1 and g2`, and `g1`. Gene IDs are normalized to lower case and `and/or`
#' are interpreted as logical separators. More complex nested rules are reserved
#' for later versions.
#'
#' @param gpr A single GPR string.
#'
#' @return A list of AND groups. The top-level list represents OR alternatives;
#' each character vector contains genes in one AND complex.
#' @export
rc_parse_gpr_simple <- function(gpr) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) {
    return(list())
  }
  gpr <- gsub("\\s+", " ", trimws(tolower(gpr)))
  gpr <- gsub("\\(", "", gpr)
  gpr <- gsub("\\)", "", gpr)

  or_parts <- strsplit(gpr, "\\s+or\\s+")[[1]]
  lapply(or_parts, function(x) {
    genes <- strsplit(x, "\\s+and\\s+")[[1]]
    genes <- trimws(genes)
    genes[nzchar(genes)]
  })
}

#' Parse a reaction GPR table
#'
#' @param gpr_table A data.frame with columns `reaction_id` and `gpr`.
#'
#' @return A named list of parsed GPR rules, keyed by reaction ID.
#' @export
rc_parse_gpr_table <- function(gpr_table) {
  if (!is.data.frame(gpr_table)) {
    stop("`gpr_table` must be a data.frame.", call. = FALSE)
  }
  missing_cols <- setdiff(c("reaction_id", "gpr"), colnames(gpr_table))
  if (length(missing_cols) > 0) {
    stop("`gpr_table` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (anyNA(gpr_table$reaction_id) || any(!nzchar(trimws(gpr_table$reaction_id)))) {
    stop("`reaction_id` values must be non-empty.", call. = FALSE)
  }
  parsed <- lapply(gpr_table$gpr, rc_parse_gpr_simple)
  names(parsed) <- as.character(gpr_table$reaction_id)
  parsed
}

#' Compute gene promiscuity weights from parsed GPR rules
#'
#' @param gpr_list A named list of parsed GPR rules.
#' @param mode Weighting mode: `"sqrt"` for `1/sqrt(N_rxn)`, `"linear"` for
#' `1/N_rxn`, or `"none"` for no promiscuity correction.
#'
#' @return A named numeric vector of gene weights.
#' @export
rc_promiscuity_weight <- function(gpr_list, mode = c("sqrt", "linear", "none")) {
  mode <- match.arg(mode)
  genes <- unique(unlist(gpr_list, use.names = FALSE))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }

  reaction_genes <- lapply(gpr_list, function(rule) unique(unlist(rule, use.names = FALSE)))
  gene_rxn <- table(unlist(reaction_genes, use.names = FALSE))
  n_rxn <- as.numeric(gene_rxn)

  w <- switch(
    mode,
    none = rep(1, length(n_rxn)),
    sqrt = 1 / sqrt(n_rxn),
    linear = 1 / n_rxn
  )
  names(w) <- names(gene_rxn)
  w
}
