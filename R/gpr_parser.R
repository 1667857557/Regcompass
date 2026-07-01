#' Parse a simple GPR rule
#'
#' Supports flat GPR forms such as `(g1 and g2) or g3`, `g1 or g2`,
#' `g1 and g2`, and `g1`. Complex nested mixtures such as
#' `g1 and (g2 or g3)` stop rather than being silently mis-parsed.
#' @export
rc_parse_gpr_simple <- function(gpr) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) return(list())
  gpr <- gsub("\\s+", " ", trimws(tolower(gpr)))
  if (grepl("\\(.*\\(", gpr) || grepl("\\).*\\)", gpr)) {
    # Multiple flat parenthesized OR terms are allowed; true nesting is not.
    flat <- gsub("\\([^()]+\\)", "X", gpr)
    if (grepl("[()]", flat)) stop("Complex or nested GPR formulas are not supported; provide a long-table GPR.", call. = FALSE)
  }
  if (grepl("and\\s*\\(", gpr) || grepl("\\)\\s*and", gpr)) {
    stop("Complex nested AND/OR GPR formulas are not supported; provide a long-table GPR.", call. = FALSE)
  }
  gpr <- gsub("[()]", "", gpr)
  or_parts <- strsplit(gpr, "\\s+or\\s+")[[1]]
  lapply(or_parts, function(x) {
    genes <- strsplit(x, "\\s+and\\s+")[[1]]
    genes <- trimws(genes)
    genes[nzchar(genes)]
  })
}

#' Parse a reaction GPR table
#'
#' Accepts either `reaction_id + gpr` or a pre-parsed long table with
#' `reaction_id + and_group_id + gene`.
#' @export
rc_parse_gpr_table <- function(gpr_table) {
  if (!is.data.frame(gpr_table)) stop("`gpr_table` must be a data.frame.", call. = FALSE)
  if (all(c("reaction_id", "and_group_id", "gene") %in% colnames(gpr_table))) {
    if (anyNA(gpr_table$reaction_id) || anyNA(gpr_table$and_group_id) || anyNA(gpr_table$gene)) stop("Long-table GPR columns must not contain NA.", call. = FALSE)
    split_rxn <- split(gpr_table, as.character(gpr_table$reaction_id))
    return(lapply(split_rxn, function(df) {
      split_gene <- split(tolower(as.character(df$gene)), df$and_group_id)
      unname(lapply(split_gene, function(x) unique(x[nzchar(x)])))
    }))
  }
  missing_cols <- setdiff(c("reaction_id", "gpr"), colnames(gpr_table))
  if (length(missing_cols) > 0) stop("`gpr_table` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (anyNA(gpr_table$reaction_id) || any(!nzchar(trimws(gpr_table$reaction_id)))) stop("`reaction_id` values must be non-empty.", call. = FALSE)
  parsed <- lapply(gpr_table$gpr, rc_parse_gpr_simple)
  names(parsed) <- as.character(gpr_table$reaction_id)
  parsed
}

#' Compute gene promiscuity weights from parsed GPR rules
#' @export
rc_promiscuity_weight <- function(gpr_list, mode = c("sqrt", "linear", "none")) {
  mode <- match.arg(mode)
  genes <- unique(unlist(gpr_list, use.names = FALSE))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) return(stats::setNames(numeric(0), character(0)))
  reaction_genes <- lapply(gpr_list, function(rule) unique(unlist(rule, use.names = FALSE)))
  gene_rxn <- table(unlist(reaction_genes, use.names = FALSE))
  n_rxn <- as.numeric(gene_rxn)
  w <- switch(mode, none = rep(1, length(n_rxn)), sqrt = 1 / sqrt(n_rxn), linear = 1 / n_rxn)
  names(w) <- names(gene_rxn)
  w
}
