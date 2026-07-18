#' Parse a Boolean GPR rule
#'
#' Parses nested `AND`/`OR` expressions with `AND` precedence and converts the
#' rule to disjunctive normal form. Each returned list element is one alternative
#' enzyme complex and contains the genes required by that complex.
#'
#' @param gpr One GPR expression.
#' @param max_terms Maximum number of DNF alternatives allowed during expansion.
#' @param preserve_case Retain source gene-symbol capitalization. The default
#'   lower-cases computational identifiers; Mouse-GEM import opts into native
#'   mouse symbol capitalization.
#' @return A list of character vectors.
.rc_gpr_tokenize <- function(gpr) {
  text <- gsub("([()])", " \\1 ", trimws(as.character(gpr)))
  tokens <- strsplit(gsub("\\s+", " ", text), " ", fixed = TRUE)[[1L]]
  tokens <- tokens[nzchar(tokens)]
  logical_operator <- tolower(tokens) %in% c("and", "or")
  tokens[logical_operator] <- tolower(tokens[logical_operator])
  tokens
}

.rc_gpr_ast_to_dnf <- function(node, max_terms = 10000L) {
  if (identical(node$type, "gene")) return(list(node$value))
  child_terms <- lapply(
    node$children,
    .rc_gpr_ast_to_dnf,
    max_terms = max_terms
  )
  if (identical(node$type, "or")) {
    answer <- unlist(child_terms, recursive = FALSE)
  } else {
    answer <- list(character())
    for (terms in child_terms) {
      answer <- unlist(lapply(answer, function(left) {
        lapply(terms, function(right) unique(c(left, right)))
      }), recursive = FALSE)
      if (length(answer) > max_terms) {
        stop(
          "GPR expansion exceeded `max_terms`; use a structured long-table GPR.",
          call. = FALSE
        )
      }
    }
  }
  keys <- vapply(answer, function(x) {
    paste(sort(unique(x)), collapse = "\001")
  }, character(1))
  answer[!duplicated(keys)]
}

rc_parse_gpr_simple <- function(gpr, max_terms = 10000L,
                                preserve_case = FALSE) {
  if (length(gpr) != 1L || is.na(gpr) || !nzchar(trimws(gpr))) return(list())
  if (!is.numeric(max_terms) || length(max_terms) != 1L ||
      !is.finite(max_terms) || max_terms < 1) {
    stop("`max_terms` must be one positive finite number.", call. = FALSE)
  }
  tokens <- .rc_gpr_tokenize(gpr)
  position <- 1L
  current <- function() {
    if (position <= length(tokens)) tokens[[position]] else NA_character_
  }
  consume <- function(expected = NULL) {
    token <- current()
    if (!is.null(expected) && !identical(token, expected)) {
      found <- if (is.na(token)) "<end>" else token
      stop(
        "Malformed GPR: expected `", expected, "` but found `", found, "`.",
        call. = FALSE
      )
    }
    position <<- position + 1L
    token
  }
  parse_primary <- NULL
  parse_and <- NULL
  parse_or <- NULL
  parse_primary <- function() {
    token <- current()
    if (is.na(token)) {
      stop("Malformed GPR: unexpected end of rule.", call. = FALSE)
    }
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
    stop(
      "Malformed GPR: unexpected trailing token `", current(), "`.",
      call. = FALSE
    )
  }
  lapply(.rc_gpr_ast_to_dnf(ast, as.integer(max_terms)), function(x) {
    x <- unique(x[nzchar(x)])
    if (isTRUE(preserve_case)) x else tolower(x)
  })
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
