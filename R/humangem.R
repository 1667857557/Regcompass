#' Download Human-GEM and prepare RegCompass GPR tables
#'
#' Downloads the official SysBioChalmers/Human-GEM GitHub repository archive,
#' extracts model annotations, and returns GPR tables ready for RegCompass. By
#' default genes are converted from Human-GEM Ensembl IDs to gene symbols using
#' `model/genes.tsv`.
#' @export
rc_download_humangem_gpr_table <- function(destdir = tempfile("Human-GEM-"),
                                           ref = "main",
                                           gene_format = c("symbol", "ensembl"),
                                           repo_url = "https://github.com/SysBioChalmers/Human-GEM",
                                           overwrite = FALSE,
                                           quiet = FALSE) {
  gene_format <- match.arg(gene_format)
  if (dir.exists(destdir) && length(list.files(destdir, all.files = TRUE, no.. = TRUE)) > 0L && !overwrite) {
    stop("`destdir` already exists and is not empty. Use `overwrite = TRUE` or choose another directory.", call. = FALSE)
  }
  if (dir.exists(destdir) && overwrite) unlink(destdir, recursive = TRUE, force = TRUE)
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)

  archive <- file.path(destdir, paste0("Human-GEM-", ref, ".zip"))
  archive_url <- paste0(sub("/$", "", repo_url), "/archive/refs/heads/", ref, ".zip")
  utils::download.file(archive_url, archive, mode = "wb", quiet = quiet)
  utils::unzip(archive, exdir = destdir)

  roots <- list.dirs(destdir, recursive = FALSE, full.names = TRUE)
  repo_dir <- roots[grepl("Human-GEM", basename(roots), ignore.case = TRUE)][1]
  if (is.na(repo_dir)) stop("Downloaded archive did not contain a Human-GEM repository directory.", call. = FALSE)

  out <- rc_prepare_humangem_gpr_table(repo_dir, gene_format = gene_format)
  attr(out, "repo_dir") <- repo_dir
  attr(out, "archive") <- archive
  attr(out, "source_url") <- repo_url
  attr(out, "ref") <- ref
  out
}

rc_prepare_humangem_gpr_table <- function(repo_dir, gene_format = c("symbol", "ensembl")) {
  gene_format <- match.arg(gene_format)
  model_dir <- file.path(repo_dir, "model")
  genes_tsv <- file.path(model_dir, "genes.tsv")
  reactions_tsv <- file.path(model_dir, "reactions.tsv")
  model_yml <- file.path(model_dir, "Human-GEM.yml")
  missing <- c(genes_tsv, reactions_tsv, model_yml)[!file.exists(c(genes_tsv, reactions_tsv, model_yml))]
  if (length(missing) > 0L) stop("Missing Human-GEM model files: ", paste(basename(missing), collapse = ", "), call. = FALSE)

  genes <- utils::read.delim(genes_tsv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("genes", "geneSymbols") %in% colnames(genes))) stop("`genes.tsv` must contain `genes` and `geneSymbols` columns.", call. = FALSE)
  gene_map <- stats::setNames(as.character(genes$geneSymbols), as.character(genes$genes))
  gene_map[!nzchar(gene_map) | is.na(gene_map)] <- names(gene_map)[!nzchar(gene_map) | is.na(gene_map)]

  reactions <- utils::read.delim(reactions_tsv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"rxns" %in% colnames(reactions)) stop("`reactions.tsv` must contain an `rxns` column.", call. = FALSE)
  reaction_rules <- rc_read_humangem_yml_rules(model_yml)
  reaction_rules <- reaction_rules[reaction_rules$reaction_id %in% reactions$rxns & nzchar(reaction_rules$gpr), , drop = FALSE]

  rows <- lapply(seq_len(nrow(reaction_rules)), function(i) {
    rid <- reaction_rules$reaction_id[[i]]
    rule <- reaction_rules$gpr[[i]]
    if (gene_format == "symbol") rule <- rc_replace_humangem_gene_ids(rule, gene_map)
    parsed <- tryCatch(rc_parse_gpr_simple(rule), error = function(e) list())
    if (length(parsed) == 0L) return(NULL)
    do.call(rbind, lapply(seq_along(parsed), function(j) {
      data.frame(reaction_id = rid, and_group_id = j, gene = parsed[[j]], stringsAsFactors = FALSE)
    }))
  })
  gpr_table <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(gpr_table)) gpr_table <- data.frame(reaction_id = character(), and_group_id = integer(), gene = character())
  gpr_table$gene <- toupper(gpr_table$gene)
  gpr_table <- unique(gpr_table)
  rownames(gpr_table) <- NULL

  list(
    gpr_table = gpr_table,
    metabolic_genes = sort(unique(gpr_table$gene)),
    reaction_rules = reaction_rules,
    genes = genes,
    reactions = reactions
  )
}

rc_read_humangem_yml_rules <- function(model_yml) {
  x <- readLines(model_yml, warn = FALSE)
  ids <- sub("^\\s*-\\s*id:\\s*['\"]?([^'\"]+)['\"]?.*$", "\\1", x[grepl("^\\s*-\\s*id:\\s*", x)])
  id_lines <- which(grepl("^\\s*-\\s*id:\\s*", x))
  rule_lines <- which(grepl("^\\s*gene_reaction_rule:\\s*", x))
  rxn_ids <- character(0)
  rules <- character(0)
  for (line in rule_lines) {
    id_idx <- max(which(id_lines < line), na.rm = TRUE)
    if (!is.finite(id_idx)) next
    rule <- sub("^\\s*gene_reaction_rule:\\s*", "", x[[line]])
    rule <- gsub("^['\"]|['\"]$", "", trimws(rule))
    rxn_ids <- c(rxn_ids, ids[[id_idx]])
    rules <- c(rules, rule)
  }
  data.frame(reaction_id = rxn_ids, gpr = rules, stringsAsFactors = FALSE)
}

rc_replace_humangem_gene_ids <- function(rule, gene_map) {
  ids <- names(gene_map)
  ids <- ids[order(nchar(ids), decreasing = TRUE)]
  for (id in ids) {
    rule <- gsub(paste0("\\b", id, "\\b"), gene_map[[id]], rule, perl = TRUE)
  }
  rule
}
