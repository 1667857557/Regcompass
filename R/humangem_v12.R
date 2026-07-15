#' Prepare a Human-GEM object with v1.2 reaction annotations
#'
#' This is the v1.2 preparation path. It retains Human-GEM subsystem labels and
#' KEGG, Reactome, Rhea and master-Rhea reaction identifiers required by the
#' GRN-driven meta-module expansion.
#' @export
rc_prepare_human2_gem_v12 <- function(version = "2.0.0",
                                      cache_dir = tools::R_user_dir("RegCompassR", "cache"),
                                      save_rds = file.path(cache_dir, paste0("Human2_", version, "_regcompass_v12.rds")),
                                      force_download = FALSE,
                                      allow_latest = FALSE) {
  if (identical(version, "latest") && !isTRUE(allow_latest)) {
    stop("`version = 'latest'` requires `allow_latest = TRUE`; use a pinned Human-GEM release.", call. = FALSE)
  }
  if (isTRUE(force_download) && file.exists(save_rds)) unlink(save_rds, force = TRUE)
  if (!file.exists(save_rds)) {
    ref <- if (identical(version, "latest")) "main" else paste0("v", version)
    tmp <- tempfile("Human-GEM-v12-")
    gpr <- rc_download_humangem_gpr_table(destdir = tmp, ref = ref, ref_type = "auto",
                                          gene_format = "symbol", overwrite = TRUE, quiet = TRUE)
    repo_dir <- attr(gpr, "repo_dir")
    model_yml <- file.path(repo_dir, "model", "Human-GEM.yml")
    checksum <- tools::md5sum(model_yml)[[1L]]
    gem_new <- rc_convert_humangem_yaml_to_regcompass(model_yml, version = version,
                                                       commit = ref, checksum = checksum)
    gem_new <- rc_enrich_humangem_v12_metadata(gem_new,
                                                reactions_tsv = gpr$reactions,
                                                model_yml = model_yml)
    gem_new <- rc_annotate_reaction_roles(gem_new, overwrite_existing = TRUE)
    gem_new$gpr_table <- gpr$gpr_table
    gem_new$metabolic_genes <- gpr$metabolic_genes
    gem_new$reaction_rules <- gpr$reaction_rules
    gem_new$genes <- gpr$genes
    gem_new$reactions <- gpr$reactions
    gem_new$model_info$gene_format <- "symbol"
    gem_new$model_info$archive <- attr(gpr, "archive")
    gem_new$model_info$archive_url <- attr(gpr, "archive_url") %||% NA_character_
    gem_new$model_info$annotation_schema <- "regcompass_humangem_v12"
    rc_validate_human2_gem(gem_new)
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(gem_new, save_rds)
  }
  gem <- rc_read_gem(save_rds)
  rc_validate_human2_gem(gem)
  gem
}

#' Add Human-GEM subsystem and cross-database reaction annotations
#' @export
rc_enrich_humangem_v12_metadata <- function(gem, reactions_tsv = NULL, model_yml = NULL) {
  gem <- rc_enrich_humangem_metadata(gem, reactions_tsv = reactions_tsv)
  gv <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  meta$reaction_id <- gv$reactions

  if (!is.null(model_yml)) {
    if (!file.exists(model_yml)) stop("Human-GEM YAML file not found: ", model_yml, call. = FALSE)
    if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required to parse Human-GEM YAML metadata.", call. = FALSE)
    model <- yaml::read_yaml(model_yml)
    rxns <- model$reactions %||% list()
    yml_ids <- vapply(rxns, function(x) as.character(x$id %||% x$reaction_id %||% NA_character_), character(1))
    collapse_field <- function(x, candidates) {
      value <- NULL
      for (nm in candidates) {
        if (!is.null(x[[nm]])) {
          value <- x[[nm]]
          break
        }
      }
      if (is.null(value)) return(NA_character_)
      value <- unlist(value, recursive = TRUE, use.names = FALSE)
      value <- .rc_mm_trim_unique(value)
      if (!length(value)) NA_character_ else paste(value, collapse = ";")
    }
    yml_meta <- data.frame(
      reaction_id = yml_ids,
      subsystem_yml = vapply(rxns, collapse_field, character(1),
                             candidates = c("subsystem", "subsystems", "subSystem", "subSystems")),
      equation_yml = vapply(rxns, collapse_field, character(1),
                            candidates = c("equation", "reaction_formula", "formula")),
      stringsAsFactors = FALSE
    )
    yml_meta <- yml_meta[!is.na(yml_meta$reaction_id) & nzchar(yml_meta$reaction_id), , drop = FALSE]
    idx <- match(meta$reaction_id, yml_meta$reaction_id)
    yml_sub <- yml_meta$subsystem_yml[idx]
    if (!"subsystem" %in% colnames(meta)) meta$subsystem <- NA_character_
    replace_sub <- (is.na(meta$subsystem) | !nzchar(as.character(meta$subsystem)) |
                    as.character(meta$subsystem) == "UNASSIGNED") & !is.na(yml_sub) & nzchar(yml_sub)
    meta$subsystem[replace_sub] <- yml_sub[replace_sub]
    yml_eq <- yml_meta$equation_yml[idx]
    if (!"equation" %in% colnames(meta)) meta$equation <- NA_character_
    replace_eq <- (is.na(meta$equation) | !nzchar(as.character(meta$equation))) & !is.na(yml_eq) & nzchar(yml_eq)
    meta$equation[replace_eq] <- yml_eq[replace_eq]
  }

  if (!is.null(reactions_tsv) && is.data.frame(reactions_tsv)) {
    id_col <- .rc_mm_first_column(reactions_tsv, c("rxns", "reaction_id", "id"))
    if (!is.null(id_col)) {
      tab <- reactions_tsv[match(meta$reaction_id, as.character(reactions_tsv[[id_col]])), , drop = FALSE]
      copy_annotation <- function(target, candidates) {
        src <- .rc_mm_first_column(tab, candidates)
        if (!is.null(src)) meta[[target]] <<- as.character(tab[[src]])
      }
      copy_annotation("kegg_reaction_id", c("rxnKEGGID", "kegg_reaction_id", "kegg_id"))
      copy_annotation("reactome_reaction_id", c("rxnREACTOMEID", "reactome_reaction_id", "reactome_id"))
      copy_annotation("rhea_reaction_id", c("rxnRheaID", "rhea_reaction_id", "rhea_id"))
      copy_annotation("rhea_master_id", c("rxnRheaMasterID", "rhea_master_id", "master_rhea_id"))
    }
  }

  if (!"subsystem" %in% colnames(meta)) meta$subsystem <- NA_character_
  meta$subsystem <- as.character(meta$subsystem)
  meta$metabolic_module <- meta$subsystem
  meta$metabolic_module[is.na(meta$metabolic_module) | !nzchar(meta$metabolic_module)] <- "UNASSIGNED"
  gem$reaction_meta <- meta
  gem
}

#' Build normalized reaction annotation maps used by meta-module expansion
#' @export
rc_reaction_crossref_maps <- function(gem, subsystem_table = NULL) {
  gv <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta)) meta <- data.frame(reaction_id = gv$reactions, stringsAsFactors = FALSE)
  if (!"reaction_id" %in% colnames(meta)) meta$reaction_id <- gv$reactions
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  meta$reaction_id <- gv$reactions

  if (!is.null(subsystem_table)) {
    if (!is.data.frame(subsystem_table)) stop("`subsystem_table` must be a data.frame.", call. = FALSE)
    rid <- .rc_mm_first_column(subsystem_table, c("reaction_id", "rxns", "id"))
    sid <- .rc_mm_first_column(subsystem_table, c("subsystem", "subSystems", "subsystem_id", "metabolic_module"))
    if (is.null(rid) || is.null(sid)) stop("`subsystem_table` must contain reaction and subsystem columns.", call. = FALSE)
    external <- subsystem_table[, c(rid, sid), drop = FALSE]
    colnames(external) <- c("reaction_id", "subsystem")
  } else {
    sid <- .rc_mm_first_column(meta, c("subsystem", "subSystems", "metabolic_module"))
    external <- if (is.null(sid)) data.frame(reaction_id = character(), subsystem = character()) else
      data.frame(reaction_id = meta$reaction_id, subsystem = as.character(meta[[sid]]), stringsAsFactors = FALSE)
  }

  expand_map <- function(ids, values, value_name, drop_unassigned = FALSE) {
    rows <- lapply(seq_along(ids), function(i) {
      vals <- .rc_mm_split_values(values[[i]])
      if (drop_unassigned) vals <- vals[!toupper(vals) %in% c("UNASSIGNED", "NA", "NONE")]
      if (!length(vals)) return(NULL)
      data.frame(reaction_id = rep(ids[[i]], length(vals)), value = vals, stringsAsFactors = FALSE)
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    out <- if (length(rows)) do.call(rbind, rows) else data.frame(reaction_id = character(), value = character(), stringsAsFactors = FALSE)
    colnames(out)[[2L]] <- value_name
    unique(out)
  }

  find_values <- function(candidates) {
    col <- .rc_mm_first_column(meta, candidates)
    if (is.null(col)) rep(NA_character_, nrow(meta)) else as.character(meta[[col]])
  }

  list(
    subsystem = expand_map(as.character(external$reaction_id), as.character(external$subsystem), "subsystem_id", TRUE),
    kegg = expand_map(meta$reaction_id, find_values(c("kegg_reaction_id", "rxnKEGGID", "kegg_id")), "kegg_id"),
    reactome = expand_map(meta$reaction_id, find_values(c("reactome_reaction_id", "rxnREACTOMEID", "reactome_id")), "reactome_id"),
    rhea = expand_map(meta$reaction_id, find_values(c("rhea_reaction_id", "rxnRheaID", "rhea_id")), "rhea_id"),
    rhea_master = expand_map(meta$reaction_id, find_values(c("rhea_master_id", "rxnRheaMasterID", "master_rhea_id")), "rhea_master_id")
  )
}
