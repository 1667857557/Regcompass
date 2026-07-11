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
  if (isTRUE(force_download) && file.exists(save_rds)) unlink(save_rds, force = TRUE)
  if (!file.exists(save_rds)) {
    ref <- if (identical(version, "latest")) "main" else paste0("v", version)
    tmp <- tempfile("Human-GEM-")
    gpr <- rc_download_humangem_gpr_table(destdir = tmp, ref = ref, ref_type = "auto", gene_format = "symbol", overwrite = TRUE, quiet = TRUE)
    repo_dir <- attr(gpr, "repo_dir")
    model_yml <- file.path(repo_dir, "model", "Human-GEM.yml")
    checksum <- tools::md5sum(model_yml)[[1]]
    gem_new <- rc_convert_humangem_yaml_to_regcompass(model_yml, version = version, commit = ref, checksum = checksum)
    gem_new <- rc_enrich_humangem_metadata(gem_new, reactions_tsv = gpr$reactions)
    gem_new <- rc_annotate_reaction_roles(gem_new, overwrite_existing = TRUE)
    gem_new$gpr_table <- gpr$gpr_table
    gem_new$metabolic_genes <- gpr$metabolic_genes
    gem_new$reaction_rules <- gpr$reaction_rules
    gem_new$genes <- gpr$genes
    gem_new$reactions <- gpr$reactions
    gem_new$model_info$gene_format <- "symbol"
    gem_new$model_info$archive <- attr(gpr, "archive")
    gem_new$model_info$archive_url <- attr(gpr, "archive_url") %||% NA_character_
    rc_validate_human2_gem(gem_new)
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
                                           ref_type = c("auto", "heads", "tags"),
                                           gene_format = c("symbol", "ensembl"),
                                           repo_url = "https://github.com/SysBioChalmers/Human-GEM",
                                           overwrite = FALSE,
                                           quiet = FALSE,
                                           download_fun = utils::download.file) {
  ref_type <- match.arg(ref_type)
  gene_format <- match.arg(gene_format)
  if (dir.exists(destdir) && length(list.files(destdir, all.files = TRUE, no.. = TRUE)) > 0L && !overwrite) {
    stop("`destdir` already exists and is not empty. Use `overwrite = TRUE` or choose another directory.", call. = FALSE)
  }
  if (dir.exists(destdir) && overwrite) unlink(destdir, recursive = TRUE, force = TRUE)
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)

  archive <- file.path(destdir, paste0("Human-GEM-", ref, ".zip"))
  part <- paste0(archive, ".part")
  make_url <- function(type) paste0(sub("/$", "", repo_url), "/archive/refs/", type, "/", ref, ".zip")
  is_semver <- grepl("^v?[0-9]+\\.[0-9]+\\.[0-9]+([.-].*)?$", ref)
  candidate_types <- switch(ref_type,
                            heads = "heads",
                            tags = "tags",
                            auto = if (is_semver) c("tags", "heads") else c("heads", "tags"))
  urls <- vapply(candidate_types, make_url, character(1))

  validate_archive <- function(path) {
    if (!file.exists(path)) return("file_missing")
    info <- file.info(path)
    if (is.na(info$size) || info$size <= 0) return("empty_file")
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    magic <- readBin(con, what = "raw", n = 4L)
    if (length(magic) < 4L || !identical(magic, charToRaw("PK\003\004"))) return("invalid_zip_magic")
    listing <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) e)
    if (inherits(listing, "error")) return(paste0("unzip_list_failed: ", conditionMessage(listing)))
    files <- gsub("^([^/]+/)", "", as.character(listing$Name))
    required <- c("model/Human-GEM.yml", "model/genes.tsv", "model/reactions.tsv")
    missing <- setdiff(required, files)
    if (length(missing)) return(paste0("missing_model_files: ", paste(missing, collapse = ", ")))
    "ok"
  }

  attempts <- vector("list", length(urls))
  archive_url <- NA_character_
  for (i in seq_along(urls)) {
    u <- urls[[i]]
    if (file.exists(part)) unlink(part, force = TRUE)
    warnings <- character()
    err <- NA_character_
    status_code <- NA_integer_
    status_code <- tryCatch(
      withCallingHandlers(
        download_fun(u, part, mode = "wb", quiet = quiet),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        err <<- conditionMessage(e)
        NA_integer_
      }
    )
    if (is.null(status_code)) status_code <- 0L
    status_ok <- identical(as.integer(status_code), 0L)
    validation <- if (status_ok) validate_archive(part) else "download_status_nonzero"
    attempts[[i]] <- data.frame(url = u,
                                download_status = as.character(status_code),
                                warning = paste(warnings, collapse = " | "),
                                error = err,
                                archive_validation = validation,
                                stringsAsFactors = FALSE)
    if (status_ok && identical(validation, "ok")) {
      if (file.exists(archive)) unlink(archive, force = TRUE)
      ok <- file.rename(part, archive)
      if (!ok) {
        file.copy(part, archive, overwrite = TRUE)
        unlink(part, force = TRUE)
      }
      archive_url <- u
      break
    }
  }
  if (file.exists(part)) unlink(part, force = TRUE)
  diagnostics <- do.call(rbind, attempts[!vapply(attempts, is.null, logical(1))])
  if (is.na(archive_url)) {
    lines <- apply(diagnostics, 1L, function(x) paste0("URL=", x[["url"]], "; download_status=", x[["download_status"]], "; warning=", x[["warning"]], "; error=", x[["error"]], "; archive_validation=", x[["archive_validation"]]))
    stop("Failed to download a valid Human-GEM archive for ref: ", ref, "\n", paste(lines, collapse = "\n"), call. = FALSE)
  }
  utils::unzip(archive, exdir = destdir)

  roots <- list.dirs(destdir, recursive = FALSE, full.names = TRUE)
  repo_dir <- roots[grepl("Human-GEM", basename(roots), ignore.case = TRUE)][1]
  if (is.na(repo_dir)) stop("Downloaded archive did not contain a Human-GEM repository directory.", call. = FALSE)

  out <- rc_prepare_humangem_gpr_table(repo_dir, gene_format = gene_format)
  attr(out, "repo_dir") <- repo_dir
  attr(out, "archive") <- archive
  attr(out, "source_url") <- repo_url
  attr(out, "ref") <- ref
  attr(out, "ref_type") <- ref_type
  attr(out, "archive_url") <- archive_url
  attr(out, "download_diagnostics") <- diagnostics
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
