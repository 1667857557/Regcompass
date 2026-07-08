#' Download Human-GEM and prepare RegCompass GPR tables
#'
#' Downloads the official SysBioChalmers/Human-GEM GitHub repository archive,
#' extracts model annotations, and returns GPR tables ready for RegCompass. By
#' default genes are converted from Human-GEM Ensembl IDs to gene symbols using
#' `model/genes.tsv`.
#' @export
rc_prepare_human2_gem <- function(version = "2.0.0",
                                  cache_dir = tools::R_user_dir("RegCompassR", "cache"),
                                  save_rds = file.path(cache_dir, paste0("Human2_", version, "_regcompass.rds")),
                                  force_download = FALSE,
                                  allow_latest = FALSE,
                                  require_model_info = TRUE) {
  if (identical(version, "latest") && !isTRUE(allow_latest)) {
    stop("`version = 'latest'` requires `allow_latest = TRUE`; use a pinned Human2 release.", call. = FALSE)
  }
  if (!file.exists(save_rds)) {
    ref <- if (identical(version, "latest")) "main" else paste0("v", version)
    tmp <- tempfile("Human-GEM-")
    gpr <- rc_download_humangem_gpr_table(destdir = tmp, ref = ref, overwrite = TRUE, quiet = TRUE)
    repo_dir <- attr(gpr, "repo_dir")
    model_yml <- file.path(repo_dir, "model", "Human-GEM.yml")
    checksum <- tools::md5sum(model_yml)[[1]]
    gem_new <- rc_convert_humangem_yaml_to_regcompass(model_yml, version = version, commit = ref, checksum = checksum)
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(gem_new, save_rds)
  }
  gem <- rc_read_gem(save_rds, require_model_info = require_model_info)
  rc_validate_human2_gem(gem)
  gem
}

#' @export
rc_validate_human2_gem <- function(gem) {
  rc_validate_gem(gem)
  if (is.null(gem$model_info)) stop("Human2 GEM requires `model_info`.", call. = FALSE)
  req_info <- c("source", "version", "commit", "checksum", "conversion_date")
  miss_info <- setdiff(req_info, names(gem$model_info))
  if (length(miss_info)) stop("Human2 `model_info` missing: ", paste(miss_info, collapse = ", "), call. = FALSE)
  if (is.null(gem$gpr_table)) stop("Human2 GEM requires `gpr_table`.", call. = FALSE)
  req_gpr <- c("reaction_id", "and_group_id", "gene")
  miss_gpr <- setdiff(req_gpr, colnames(gem$gpr_table))
  if (length(miss_gpr)) stop("Human2 `gpr_table` missing columns: ", paste(miss_gpr, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

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
  rule_pattern <- "^\\s*-?\\s*gene_reaction_rule:\\s*"
  rule_lines <- which(grepl(rule_pattern, x))
  rxn_ids <- character(0)
  rules <- character(0)
  for (line in rule_lines) {
    id_idx <- max(which(id_lines < line), na.rm = TRUE)
    if (!is.finite(id_idx)) next
    rule <- sub(rule_pattern, "", x[[line]])
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

#' Convert a Human-GEM YAML model to a RegCompass GEM object
#' @export
rc_convert_humangem_yaml_to_regcompass <- function(model_yml,
                                                   version = NA_character_,
                                                   commit = NA_character_,
                                                   checksum = NA_character_) {
  if (!file.exists(model_yml)) stop("Human-GEM YAML file not found: ", model_yml, call. = FALSE)
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required to parse Human-GEM YAML files.", call. = FALSE)
  model <- yaml::read_yaml(model_yml)
  mets <- model$metabolites %||% list()
  rxns <- model$reactions %||% list()
  metabolite_ids <- vapply(mets, function(x) as.character(x$id %||% x$metabolite_id %||% NA_character_), character(1))
  reaction_ids <- vapply(rxns, function(x) as.character(x$id %||% x$reaction_id %||% NA_character_), character(1))
  metabolite_ids <- metabolite_ids[!is.na(metabolite_ids) & nzchar(metabolite_ids)]
  reaction_ids <- reaction_ids[!is.na(reaction_ids) & nzchar(reaction_ids)]
  S <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(length(metabolite_ids), length(reaction_ids)), dimnames = list(metabolite_ids, reaction_ids))
  lb <- stats::setNames(rep(-1000, length(reaction_ids)), reaction_ids)
  ub <- stats::setNames(rep(1000, length(reaction_ids)), reaction_ids)
  gpr_rows <- list()
  for (j in seq_along(rxns)) {
    rid <- reaction_ids[[j]]
    r <- rxns[[j]]
    if (!is.null(r$lower_bound)) lb[[rid]] <- as.numeric(r$lower_bound)
    if (!is.null(r$upper_bound)) ub[[rid]] <- as.numeric(r$upper_bound)
    sto <- r$metabolites %||% r$stoichiometry
    if (length(sto)) {
      ids <- names(sto); vals <- as.numeric(unlist(sto, use.names = FALSE)); idx <- match(ids, metabolite_ids); ok <- !is.na(idx) & is.finite(vals)
      if (any(ok)) S[cbind(idx[ok], j)] <- vals[ok]
    }
    rule <- as.character(r$gene_reaction_rule %||% r$grRule %||% "")
    if (nzchar(rule)) {
      parsed <- tryCatch(rc_parse_gpr_simple(rule), error = function(e) list())
      if (length(parsed)) gpr_rows[[length(gpr_rows) + 1L]] <- do.call(rbind, lapply(seq_along(parsed), function(k) data.frame(reaction_id = rid, and_group_id = k, gene = parsed[[k]], stringsAsFactors = FALSE)))
    }
  }
  reaction_meta <- data.frame(reaction_id = reaction_ids, name = vapply(rxns, function(x) as.character(x$name %||% NA_character_), character(1)), stringsAsFactors = FALSE)
  metabolite_meta <- data.frame(metabolite_id = metabolite_ids, name = vapply(mets, function(x) as.character(x$name %||% NA_character_), character(1)), stringsAsFactors = FALSE)
  gpr_table <- if (length(gpr_rows)) unique(do.call(rbind, gpr_rows)) else data.frame(reaction_id = character(), and_group_id = integer(), gene = character())
  out <- list(S = S, lb = lb, ub = ub, reaction_meta = reaction_meta, metabolite_meta = metabolite_meta, gpr_table = gpr_table, model_info = list(source = "SysBioChalmers/Human-GEM", version = version, commit = commit, checksum = checksum, conversion_date = as.character(Sys.Date())))
  rc_validate_gem(out)
  out
}
