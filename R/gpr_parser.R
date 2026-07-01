#' Parse a simple GPR rule
#'
#' Supports `g1`, `g1 and g2`, `g1 or g2`, and flat OR-of-AND rules such as
#' `(g1 and g2) or (g3 and g4)`. For complex nested rules, provide a long table
#' with `reaction_id`, `and_group_id`, and `gene`.
#' @export
rc_parse_gpr_simple <- function(gpr) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) return(list())
  gpr <- gsub("\\s+", " ", trimws(tolower(gpr)))
  if (grepl("and\\s*\\([^)]*\\s+or\\s+", gpr) || grepl("or\\s*\\([^)]*\\s+and\\s*\\([^)]*", gpr)) {
    stop("Complex nested GPR formulas are not supported; provide a long-table GPR.", call. = FALSE)
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
  if (length(missing_cols) > 0L) stop("`gpr_table` is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (anyNA(gpr_table$reaction_id) || any(!nzchar(trimws(gpr_table$reaction_id)))) stop("`reaction_id` values must be non-empty.", call. = FALSE)
  parsed <- lapply(gpr_table$gpr, rc_parse_gpr_simple)
  names(parsed) <- as.character(gpr_table$reaction_id)
  parsed
}

rc_promiscuity_weight <- function(gpr_list) {
  genes <- unique(unlist(gpr_list, use.names = FALSE))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0L) return(stats::setNames(numeric(0), character(0)))
  reaction_genes <- lapply(gpr_list, function(rule) unique(unlist(rule, use.names = FALSE)))
  tab <- table(unlist(reaction_genes, use.names = FALSE))
  n_rxn <- as.numeric(tab)
  w <- 1 / sqrt(n_rxn)
  names(w) <- names(tab)
  w
}
