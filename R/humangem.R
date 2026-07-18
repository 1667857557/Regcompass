.rc_species_gem_spec <- function(species = c("human", "mouse"), version = NULL) {
  species <- match.arg(species)
  if (is.null(version) || length(version) != 1L || is.na(version) ||
      !nzchar(trimws(as.character(version)))) {
    version <- if (identical(species, "human")) "2.0.0" else "1.8.0"
  }
  if (!is.character(version) || length(version) != 1L ||
      is.na(version) || !nzchar(trimws(version))) {
    stop("`version` must be one non-empty character value.", call. = FALSE)
  }
  version <- trimws(version)
  if (identical(species, "human")) {
    return(list(
      species = "human",
      species_label = "Homo sapiens",
      taxonomy_id = "9606",
      source = "SysBioChalmers/Human-GEM",
      repo_url = "https://github.com/SysBioChalmers/Human-GEM",
      repository_name = "Human-GEM",
      model_file = "Human-GEM.yml",
      version = version,
      default_version = "2.0.0",
      cache_prefix = "Human2",
      gene_format = "symbol",
      citation = "Robinson et al., Sci Signal 2020; Human-GEM consortium",
      citation_doi = "10.1126/scisignal.aaz1482"
    ))
  }
  list(
    species = "mouse",
    species_label = "Mus musculus",
    taxonomy_id = "10090",
    source = "SysBioChalmers/Mouse-GEM",
    repo_url = "https://github.com/SysBioChalmers/Mouse-GEM",
    repository_name = "Mouse-GEM",
    model_file = "Mouse-GEM.yml",
    version = version,
    default_version = "1.8.0",
    cache_prefix = "Mouse",
    gene_format = "symbol",
    citation = "Wang et al., PNAS 2021",
    citation_doi = "10.1073/pnas.2102344118"
  )
}

.rc_load_compatible_species_gem <- function(save_rds, spec) {
  if (!file.exists(save_rds)) return(NULL)
  cached <- tryCatch(
    rc_read_gem(save_rds),
    error = function(error) error
  )
  reason <- NULL
  if (inherits(cached, "error")) {
    reason <- conditionMessage(cached)
  } else {
    reason <- tryCatch(
      {
        rc_validate_species_gem(cached, spec$species)
        recorded_source <- as.character(cached$model_info$source %||% "")
        recorded_version <- as.character(cached$model_info$version %||% "")
        if (!identical(recorded_source, spec$source)) {
          stop(
            "cached model source is `", recorded_source,
            "` instead of `", spec$source, "`",
            call. = FALSE
          )
        }
        if (!identical(recorded_version, spec$version)) {
          stop(
            "cached model version is `", recorded_version,
            "` instead of `", spec$version, "`",
            call. = FALSE
          )
        }
        NULL
      },
      error = function(error) conditionMessage(error)
    )
  }
  if (!is.null(reason)) {
    warning(
      "Removing incompatible cached ", spec$repository_name,
      " model at `", save_rds, "`: ", reason,
      call. = FALSE
    )
    unlink(save_rds, force = TRUE)
    return(NULL)
  }
  cached
}

#' Prepare a species-specific genome-scale metabolic model
#'
#' Downloads and converts a pinned official SysBioChalmers GEM release. Human
#' mode uses Human-GEM and mouse mode uses Mouse-GEM directly; mouse genes are
#' retained as mouse symbols and are not converted through human orthologues.
#'
#' @param species Model organism: `"human"` or `"mouse"`.
#' @param version Pinned release version. Defaults to Human-GEM 2.0.0 or
#'   Mouse-GEM 1.8.0.
#' @param cache_dir Persistent model cache directory.
#' @param save_rds Optional output RDS path. When `NULL`, a species/version path
#'   is generated inside `cache_dir`.
#' @param force_download Re-download and rebuild an existing cached model.
#' @param allow_latest Permit the unpinned `version = "latest"` mode.
#' @return A validated RegCompass GEM list.
#' @export
rc_prepare_gem <- function(
    species = c("human", "mouse"),
    version = NULL,
    cache_dir = tools::R_user_dir("RegCompassR", "cache"),
    save_rds = NULL,
    force_download = FALSE,
    allow_latest = FALSE) {
  species <- match.arg(species)
  spec <- .rc_species_gem_spec(species, version)
  if (identical(spec$version, "latest") && !isTRUE(allow_latest)) {
    stop(
      "`version = 'latest'` requires `allow_latest = TRUE`; use a pinned GEM release for reproducible analyses.",
      call. = FALSE
    )
  }
  if (!is.logical(force_download) || length(force_download) != 1L ||
      is.na(force_download) ||
      !is.logical(allow_latest) || length(allow_latest) != 1L ||
      is.na(allow_latest)) {
    stop("`force_download` and `allow_latest` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(save_rds)) {
    save_rds <- file.path(
      cache_dir,
      paste0(spec$cache_prefix, "_", spec$version, "_regcompass.rds")
    )
  }
  if (!is.character(save_rds) || length(save_rds) != 1L ||
      is.na(save_rds) || !nzchar(save_rds)) {
    stop("`save_rds` must be one non-empty file path or NULL.", call. = FALSE)
  }
  if (isTRUE(force_download) && file.exists(save_rds)) {
    unlink(save_rds, force = TRUE)
  }
  cached <- if (isTRUE(force_download)) {
    NULL
  } else {
    .rc_load_compatible_species_gem(save_rds, spec)
  }
  if (!is.null(cached)) return(cached)

  ref <- if (identical(spec$version, "latest")) {
    "main"
  } else {
    paste0("v", spec$version)
  }
  tmp <- tempfile(paste0(spec$repository_name, "-"))
  prepared <- rc_download_species_gem(
    species = species,
    destdir = tmp,
    ref = ref,
    ref_type = "auto",
    gene_format = spec$gene_format,
    overwrite = TRUE,
    quiet = TRUE
  )
  repo_dir <- attr(prepared, "repo_dir")
  model_yml <- file.path(repo_dir, "model", spec$model_file)
  checksum <- unname(tools::md5sum(model_yml)[[1L]])
  gem <- rc_convert_yaml_to_regcompass(
    model_yml = model_yml,
    species = species,
    version = spec$version,
    commit = ref,
    checksum = checksum
  )
  gem <- rc_enrich_humangem_metadata(
    gem,
    reactions_tsv = prepared$reactions,
    model_yml = model_yml
  )
  gem <- rc_annotate_reaction_roles(gem, overwrite_existing = TRUE)
  gem$gpr_table <- prepared$gpr_table
  gem$metabolic_genes <- prepared$metabolic_genes
  gem$reaction_rules <- prepared$reaction_rules
  gem$genes <- prepared$genes
  gem$reactions <- prepared$reactions
  gem$model_info$gene_format <- spec$gene_format
  gem$model_info$archive <- attr(prepared, "archive")
  gem$model_info$archive_url <- attr(prepared, "archive_url") %||%
    NA_character_
  gem$model_info$annotation_schema <- "regcompass_species_gem_v1"
  gem$model_info$citation <- spec$citation
  gem$model_info$citation_doi <- spec$citation_doi
  rc_validate_species_gem(gem, species)
  dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gem, save_rds)
  gem <- rc_read_gem(save_rds)
  rc_validate_species_gem(gem, species)
  gem
}

#' Prepare Human-GEM 2 for RegCompass
#'
#' Species-specific Human-GEM 2 entry point. Use this named helper when the
#' analysis should explicitly download, cache, and validate the Human-GEM 2
#' model path; it delegates to `rc_prepare_gem(species = "human")` with
#' Human-GEM 2 defaults.
#'
#' @inheritParams rc_prepare_gem
#' @param version Human-GEM 2 release version.
#' @export
rc_prepare_human2_gem <- function(
    version = "2.0.0",
    cache_dir = tools::R_user_dir("RegCompassR", "cache"),
    save_rds = file.path(
      cache_dir,
      paste0("Human2_", version, "_regcompass.rds")
    ),
    force_download = FALSE,
    allow_latest = FALSE) {
  rc_prepare_gem(
    species = "human",
    version = version,
    cache_dir = cache_dir,
    save_rds = save_rds,
    force_download = force_download,
    allow_latest = allow_latest
  )
}

#' Prepare Mouse-GEM for RegCompass
#'
#' Species-specific convenience entry point for users who want the Mouse-GEM
#' download, cache naming, and validation path explicitly. Equivalent to
#' `rc_prepare_gem(species = "mouse")` with Mouse-GEM defaults.
#'
#' @inheritParams rc_prepare_gem
#' @param version Mouse-GEM release version.
#' @export
rc_prepare_mouse_gem <- function(
    version = "1.8.0",
    cache_dir = tools::R_user_dir("RegCompassR", "cache"),
    save_rds = file.path(
      cache_dir,
      paste0("Mouse_", version, "_regcompass.rds")
    ),
    force_download = FALSE,
    allow_latest = FALSE) {
  rc_prepare_gem(
    species = "mouse",
    version = version,
    cache_dir = cache_dir,
    save_rds = save_rds,
    force_download = force_download,
    allow_latest = allow_latest
  )
}

rc_validate_species_gem <- function(gem, species = c("human", "mouse")) {
  species <- match.arg(species)
  rc_validate_gem(gem)
  if (is.null(gem$model_info)) {
    stop("Species GEM requires `model_info`.", call. = FALSE)
  }
  required_info <- c(
    "source", "species", "species_label", "taxonomy_id", "version",
    "commit", "checksum", "conversion_date"
  )
  missing_info <- setdiff(required_info, names(gem$model_info))
  if (length(missing_info)) {
    stop(
      "Species GEM `model_info` missing: ",
      paste(missing_info, collapse = ", "),
      call. = FALSE
    )
  }
  if (!identical(as.character(gem$model_info$species), species)) {
    stop(
      "GEM species mismatch: expected `", species, "` but model records `",
      as.character(gem$model_info$species), "`.",
      call. = FALSE
    )
  }
  if (is.null(gem$gpr_table)) {
    stop("Species GEM requires `gpr_table`.", call. = FALSE)
  }
  required_gpr <- c("reaction_id", "and_group_id", "gene")
  missing_gpr <- setdiff(required_gpr, colnames(gem$gpr_table))
  if (length(missing_gpr)) {
    stop(
      "Species GEM `gpr_table` missing columns: ",
      paste(missing_gpr, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

rc_download_species_gem <- function(
    species = c("human", "mouse"),
    destdir = tempfile("species-GEM-"),
    ref = "main",
    ref_type = c("auto", "heads", "tags"),
    gene_format = c("symbol", "ensembl"),
    overwrite = FALSE,
    quiet = FALSE,
    download_fun = utils::download.file) {
  species <- match.arg(species)
  ref_type <- match.arg(ref_type)
  gene_format <- match.arg(gene_format)
  spec <- .rc_species_gem_spec(species, if (identical(ref, "main")) "latest" else sub("^v", "", ref))
  if (dir.exists(destdir) &&
      length(list.files(destdir, all.files = TRUE, no.. = TRUE)) > 0L &&
      !overwrite) {
    stop(
      "`destdir` already exists and is not empty. Use `overwrite = TRUE` or choose another directory.",
      call. = FALSE
    )
  }
  if (dir.exists(destdir) && overwrite) {
    unlink(destdir, recursive = TRUE, force = TRUE)
  }
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)

  archive <- file.path(destdir, paste0(spec$repository_name, "-", ref, ".zip"))
  part <- paste0(archive, ".part")
  make_url <- function(type) {
    paste0(
      sub("/$", "", spec$repo_url),
      "/archive/refs/", type, "/", ref, ".zip"
    )
  }
  is_semver <- grepl("^v?[0-9]+\\.[0-9]+\\.[0-9]+([.-].*)?$", ref)
  candidate_types <- switch(
    ref_type,
    heads = "heads",
    tags = "tags",
    auto = if (is_semver) c("tags", "heads") else c("heads", "tags")
  )
  urls <- vapply(candidate_types, make_url, character(1))
  required <- c(
    file.path("model", spec$model_file),
    file.path("model", "reactions.tsv")
  )
  if (identical(species, "human")) {
    required <- c(required, file.path("model", "genes.tsv"))
  }
  validate_archive <- function(path) {
    if (!file.exists(path)) return("file_missing")
    info <- file.info(path)
    if (is.na(info$size) || info$size <= 0) return("empty_file")
    con <- file(path, "rb")
    on.exit(close(con), add = TRUE)
    magic <- readBin(con, what = "raw", n = 4L)
    if (length(magic) < 4L ||
        !identical(magic, charToRaw("PK\003\004"))) {
      return("invalid_zip_magic")
    }
    listing <- tryCatch(utils::unzip(path, list = TRUE), error = function(e) e)
    if (inherits(listing, "error")) {
      return(paste0("unzip_list_failed: ", conditionMessage(listing)))
    }
    files <- gsub("^([^/]+/)", "", as.character(listing$Name))
    missing <- setdiff(required, files)
    if (length(missing)) {
      return(paste0("missing_model_files: ", paste(missing, collapse = ", ")))
    }
    "ok"
  }

  attempts <- vector("list", length(urls))
  archive_url <- NA_character_
  for (i in seq_along(urls)) {
    url <- urls[[i]]
    if (file.exists(part)) unlink(part, force = TRUE)
    warnings <- character()
    error_message <- NA_character_
    status_code <- tryCatch(
      withCallingHandlers(
        download_fun(url, part, mode = "wb", quiet = quiet),
        warning = function(warning) {
          warnings <<- c(warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(error) {
        error_message <<- conditionMessage(error)
        NA_integer_
      }
    )
    if (is.null(status_code)) status_code <- 0L
    status_ok <- identical(as.integer(status_code), 0L)
    validation <- if (status_ok) {
      validate_archive(part)
    } else {
      "download_status_nonzero"
    }
    attempts[[i]] <- data.frame(
      url = url,
      download_status = as.character(status_code),
      warning = paste(warnings, collapse = " | "),
      error = error_message,
      archive_validation = validation,
      stringsAsFactors = FALSE
    )
    if (status_ok && identical(validation, "ok")) {
      if (file.exists(archive)) unlink(archive, force = TRUE)
      renamed <- file.rename(part, archive)
      if (!renamed) {
        file.copy(part, archive, overwrite = TRUE)
        unlink(part, force = TRUE)
      }
      archive_url <- url
      break
    }
  }
  if (file.exists(part)) unlink(part, force = TRUE)
  diagnostics <- do.call(
    rbind,
    attempts[!vapply(attempts, is.null, logical(1))]
  )
  if (is.na(archive_url)) {
    lines <- apply(diagnostics, 1L, function(x) {
      paste0(
        "URL=", x[["url"]],
        "; download_status=", x[["download_status"]],
        "; warning=", x[["warning"]],
        "; error=", x[["error"]],
        "; archive_validation=", x[["archive_validation"]]
      )
    })
    stop(
      "Failed to download a valid ", spec$repository_name,
      " archive for ref: ", ref, "\n", paste(lines, collapse = "\n"),
      call. = FALSE
    )
  }
  utils::unzip(archive, exdir = destdir)
  roots <- list.dirs(destdir, recursive = FALSE, full.names = TRUE)
  repo_dir <- roots[
    grepl(spec$repository_name, basename(roots), ignore.case = TRUE)
  ][1L]
  if (is.na(repo_dir)) {
    stop(
      "Downloaded archive did not contain a ", spec$repository_name,
      " repository directory.",
      call. = FALSE
    )
  }
  output <- rc_prepare_species_gpr_table(
    repo_dir,
    species = species,
    gene_format = gene_format
  )
  attr(output, "repo_dir") <- repo_dir
  attr(output, "archive") <- archive
  attr(output, "source_url") <- spec$repo_url
  attr(output, "ref") <- ref
  attr(output, "ref_type") <- ref_type
  attr(output, "archive_url") <- archive_url
  attr(output, "download_diagnostics") <- diagnostics
  output
}

rc_prepare_species_gpr_table <- function(
    repo_dir,
    species = c("human", "mouse"),
    gene_format = c("symbol", "ensembl")) {
  species <- match.arg(species)
  gene_format <- match.arg(gene_format)
  spec <- .rc_species_gem_spec(species)
  model_dir <- file.path(repo_dir, "model")
  reactions_tsv <- file.path(model_dir, "reactions.tsv")
  model_yml <- file.path(model_dir, spec$model_file)
  genes_tsv <- file.path(model_dir, "genes.tsv")
  required <- c(reactions_tsv, model_yml)
  if (identical(species, "human")) required <- c(required, genes_tsv)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop(
      "Missing ", spec$repository_name, " model files: ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }
  reactions <- utils::read.delim(
    reactions_tsv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!"rxns" %in% colnames(reactions)) {
    stop("`reactions.tsv` must contain an `rxns` column.", call. = FALSE)
  }
  genes <- NULL
  gene_map <- NULL
  if (identical(species, "human")) {
    genes <- utils::read.delim(
      genes_tsv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (!all(c("genes", "geneSymbols") %in% colnames(genes))) {
      stop(
        "Human-GEM `genes.tsv` must contain `genes` and `geneSymbols` columns.",
        call. = FALSE
      )
    }
    gene_map <- stats::setNames(
      as.character(genes$geneSymbols),
      as.character(genes$genes)
    )
    bad <- is.na(gene_map) | !nzchar(gene_map)
    gene_map[bad] <- names(gene_map)[bad]
  }
  reaction_rules <- rc_read_gem_yml_rules(model_yml)
  reaction_rules <- reaction_rules[
    reaction_rules$reaction_id %in% reactions$rxns &
      nzchar(reaction_rules$gpr),
    ,
    drop = FALSE
  ]
  rows <- lapply(seq_len(nrow(reaction_rules)), function(i) {
    reaction_id <- reaction_rules$reaction_id[[i]]
    rule <- reaction_rules$gpr[[i]]
    if (identical(species, "human") && identical(gene_format, "symbol")) {
      rule <- rc_replace_gem_gene_ids(rule, gene_map)
    }
    parsed <- tryCatch(
      rc_parse_gpr_simple(
        rule, preserve_case = identical(species, "mouse")
      ),
      error = function(error) {
        stop(
          "Failed to parse ", spec$repository_name, " GPR for reaction `",
          reaction_id, "`: ", conditionMessage(error), " Rule: ", rule,
          call. = FALSE
        )
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
  gpr_table <- if (length(nonempty)) {
    do.call(rbind, nonempty)
  } else {
    data.frame(
      reaction_id = character(),
      and_group_id = integer(),
      gene = character()
    )
  }
  gpr_table$gene <- if (identical(species, "human")) {
    toupper(gpr_table$gene)
  } else {
    as.character(gpr_table$gene)
  }
  gpr_table <- unique(gpr_table)
  rownames(gpr_table) <- NULL
  if (is.null(genes)) {
    mouse_genes <- sort(unique(as.character(gpr_table$gene)))
    genes <- data.frame(
      genes = mouse_genes,
      geneSymbols = mouse_genes,
      stringsAsFactors = FALSE
    )
  }
  list(
    gpr_table = gpr_table,
    metabolic_genes = sort(unique(gpr_table$gene)),
    reaction_rules = reaction_rules,
    genes = genes,
    reactions = reactions
  )
}

rc_read_gem_yml_rules <- function(model_yml) {
  lines <- readLines(model_yml, warn = FALSE)
  id_lines <- which(grepl("^\\s*-\\s*id:\\s*", lines))
  ids <- sub(
    "^\\s*-\\s*id:\\s*['\"]?([^'\"]+)['\"]?.*$",
    "\\1",
    lines[id_lines]
  )
  rule_pattern <- "^\\s*-?\\s*gene_reaction_rule:\\s*"
  rule_lines <- which(grepl(rule_pattern, lines))
  reaction_ids <- character()
  rules <- character()
  for (line in rule_lines) {
    previous <- which(id_lines < line)
    if (!length(previous)) next
    id_index <- previous[[length(previous)]]
    rule <- sub(rule_pattern, "", lines[[line]])
    rule <- gsub("^['\"]|['\"]$", "", trimws(rule))
    reaction_ids <- c(reaction_ids, ids[[id_index]])
    rules <- c(rules, rule)
  }
  data.frame(
    reaction_id = reaction_ids,
    gpr = rules,
    stringsAsFactors = FALSE
  )
}

rc_read_humangem_yml_rules <- rc_read_gem_yml_rules

rc_replace_gem_gene_ids <- function(rule, gene_map) {
  ids <- names(gene_map)
  ids <- ids[order(nchar(ids), decreasing = TRUE)]
  for (id in ids) {
    rule <- gsub(
      paste0("\\b", id, "\\b"),
      gene_map[[id]],
      rule,
      perl = TRUE
    )
  }
  rule
}

rc_convert_yaml_to_regcompass <- function(
    model_yml,
    species = c("human", "mouse"),
    version = NA_character_,
    commit = NA_character_,
    checksum = NA_character_) {
  species <- match.arg(species)
  spec <- .rc_species_gem_spec(species, version)
  if (!file.exists(model_yml)) {
    stop(spec$repository_name, " YAML file not found: ", model_yml, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop(
      "Package 'yaml' is required to parse species GEM YAML files.",
      call. = FALSE
    )
  }
  model <- yaml::read_yaml(model_yml)
  metabolites <- model$metabolites %||% list()
  reactions <- model$reactions %||% list()
  metabolite_ids <- vapply(metabolites, function(x) {
    as.character(x$id %||% x$metabolite_id %||% NA_character_)
  }, character(1))
  reaction_ids <- vapply(reactions, function(x) {
    as.character(x$id %||% x$reaction_id %||% NA_character_)
  }, character(1))
  metabolite_ids <- metabolite_ids[
    !is.na(metabolite_ids) & nzchar(metabolite_ids)
  ]
  reaction_ids <- reaction_ids[
    !is.na(reaction_ids) & nzchar(reaction_ids)
  ]
  if (anyDuplicated(metabolite_ids) || anyDuplicated(reaction_ids)) {
    stop("Species GEM YAML contains duplicated metabolite or reaction IDs.", call. = FALSE)
  }
  S <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(length(metabolite_ids), length(reaction_ids)),
    dimnames = list(metabolite_ids, reaction_ids)
  )
  lb <- stats::setNames(rep(-1000, length(reaction_ids)), reaction_ids)
  ub <- stats::setNames(rep(1000, length(reaction_ids)), reaction_ids)
  gpr_rows <- list()
  for (j in seq_along(reactions)) {
    reaction <- reactions[[j]]
    reaction_id <- as.character(
      reaction$id %||% reaction$reaction_id %||% NA_character_
    )
    if (is.na(reaction_id) || !nzchar(reaction_id)) next
    if (!is.null(reaction$lower_bound)) {
      lb[[reaction_id]] <- as.numeric(reaction$lower_bound)
    }
    if (!is.null(reaction$upper_bound)) {
      ub[[reaction_id]] <- as.numeric(reaction$upper_bound)
    }
    stoichiometry <- reaction$metabolites %||% reaction$stoichiometry
    if (length(stoichiometry)) {
      ids <- names(stoichiometry)
      values <- as.numeric(unlist(stoichiometry, use.names = FALSE))
      row_index <- match(ids, metabolite_ids)
      valid <- !is.na(row_index) & is.finite(values)
      if (any(valid)) S[cbind(row_index[valid], j)] <- values[valid]
    }
    rule <- as.character(
      reaction$gene_reaction_rule %||% reaction$grRule %||% ""
    )
    if (nzchar(rule)) {
      parsed <- tryCatch(
        rc_parse_gpr_simple(
          rule, preserve_case = identical(species, "mouse")
        ),
        error = function(error) {
          stop(
            "Failed to parse ", spec$repository_name,
            " GPR for reaction `", reaction_id, "`: ",
            conditionMessage(error),
            call. = FALSE
          )
        }
      )
      if (length(parsed)) {
        gpr_rows[[length(gpr_rows) + 1L]] <- do.call(
          rbind,
          lapply(seq_along(parsed), function(group_index) {
            data.frame(
              reaction_id = reaction_id,
              and_group_id = group_index,
              gene = parsed[[group_index]],
              stringsAsFactors = FALSE
            )
          })
        )
      }
    }
  }
  reaction_meta <- data.frame(
    reaction_id = reaction_ids,
    name = vapply(reactions, function(x) {
      as.character(x$name %||% NA_character_)
    }, character(1)),
    stringsAsFactors = FALSE
  )
  metabolite_meta <- data.frame(
    metabolite_id = metabolite_ids,
    name = vapply(metabolites, function(x) {
      as.character(x$name %||% NA_character_)
    }, character(1)),
    compartment = vapply(metabolites, function(x) {
      as.character(x$compartment %||% NA_character_)
    }, character(1)),
    stringsAsFactors = FALSE
  )
  gpr_table <- if (length(gpr_rows)) {
    unique(do.call(rbind, gpr_rows))
  } else {
    data.frame(
      reaction_id = character(),
      and_group_id = integer(),
      gene = character()
    )
  }
  output <- list(
    S = S,
    lb = lb,
    ub = ub,
    reaction_meta = reaction_meta,
    metabolite_meta = metabolite_meta,
    gpr_table = gpr_table,
    model_info = list(
      source = spec$source,
      species = species,
      species_label = spec$species_label,
      taxonomy_id = spec$taxonomy_id,
      version = version,
      commit = commit,
      checksum = checksum,
      conversion_date = as.character(Sys.Date())
    )
  )
  rc_validate_gem(output)
  output
}
