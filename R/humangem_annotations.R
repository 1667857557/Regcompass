.rc_humangem_met_compartment_from_id <- function(metabolite_id) {
  x <- as.character(metabolite_id)
  known <- c("c", "e", "g", "l", "m", "n", "p", "r", "x")
  bracket <- sub("^.*\\[([A-Za-z])\\]$", "\\1", x)
  last <- substring(x, nchar(x), nchar(x))
  out <- ifelse(bracket %in% known, bracket, ifelse(last %in% known, last, NA_character_))
  out
}

#' Enrich Human-GEM metadata with module, equation, reversibility, and compartment annotations
.rc_enrich_humangem_metadata_core <- function(gem, reactions_tsv = NULL, metabolites_tsv = NULL) {
  gv <- rc_validate_gem(gem)

  rxn_meta <- gem$reaction_meta
  if (is.null(rxn_meta) || !is.data.frame(rxn_meta)) {
    rxn_meta <- data.frame(reaction_id = gv$reactions, stringsAsFactors = FALSE)
  }
  if (!"reaction_id" %in% colnames(rxn_meta)) rxn_meta$reaction_id <- gv$reactions
  rxn_meta <- rxn_meta[match(gv$reactions, as.character(rxn_meta$reaction_id)), , drop = FALSE]
  rxn_meta$reaction_id <- gv$reactions

  if (!is.null(reactions_tsv) && is.data.frame(reactions_tsv)) {
    id_col <- .rc_first_existing_col(reactions_tsv, c("rxns", "reaction_id", "id"))
    if (!is.null(id_col)) {
      tab <- reactions_tsv[match(gv$reactions, as.character(reactions_tsv[[id_col]])), , drop = FALSE]
      nm_col <- .rc_first_existing_col(tab, c("rxnNames", "name", "reaction_name"))
      ss_col <- .rc_first_existing_col(tab, c("subSystems", "subsystem", "sub_system"))
      eq_col <- .rc_first_existing_col(tab, c("equation", "rxnEquations", "reaction_formula"))
      rev_col <- .rc_first_existing_col(tab, c("rev", "reversible"))
      if (!is.null(nm_col)) rxn_meta$name <- as.character(tab[[nm_col]])
      if (!is.null(ss_col)) rxn_meta$subsystem <- as.character(tab[[ss_col]])
      if (!is.null(eq_col)) rxn_meta$equation <- as.character(tab[[eq_col]])
      if (!is.null(rev_col)) rxn_meta$reversible <- as.logical(tab[[rev_col]])
    }
  }
  if (!"subsystem" %in% colnames(rxn_meta)) rxn_meta$subsystem <- NA_character_
  rxn_meta$metabolic_module <- as.character(rxn_meta$subsystem)
  miss_mod <- is.na(rxn_meta$metabolic_module) | !nzchar(rxn_meta$metabolic_module)
  rxn_meta$metabolic_module[miss_mod] <- "UNASSIGNED"

  met_meta <- gem$metabolite_meta
  if (is.null(met_meta) || !is.data.frame(met_meta)) {
    met_meta <- data.frame(metabolite_id = gv$metabolites, stringsAsFactors = FALSE)
  }
  if (!"metabolite_id" %in% colnames(met_meta)) met_meta$metabolite_id <- gv$metabolites
  met_meta <- met_meta[match(gv$metabolites, as.character(met_meta$metabolite_id)), , drop = FALSE]
  met_meta$metabolite_id <- gv$metabolites

  if (!is.null(metabolites_tsv) && is.data.frame(metabolites_tsv)) {
    id_col <- .rc_first_existing_col(metabolites_tsv, c("mets", "metabolite_id", "id"))
    if (!is.null(id_col)) {
      tab <- metabolites_tsv[match(gv$metabolites, as.character(metabolites_tsv[[id_col]])), , drop = FALSE]
      comp_col <- .rc_first_existing_col(tab, c("compartment", "compartments", "metComps"))
      if (!is.null(comp_col)) met_meta$compartment <- as.character(tab[[comp_col]])
    }
  }
  if (!"compartment" %in% colnames(met_meta)) {
    met_meta$compartment <- .rc_humangem_met_compartment_from_id(met_meta$metabolite_id)
  } else {
    miss <- is.na(met_meta$compartment) | !nzchar(as.character(met_meta$compartment))
    met_meta$compartment[miss] <- .rc_humangem_met_compartment_from_id(met_meta$metabolite_id[miss])
  }

  gem$reaction_meta <- rxn_meta
  gem$metabolite_meta <- met_meta
  gem
}

rc_enrich_humangem_metadata <- function(gem, reactions_tsv = NULL, model_yml = NULL) {
  gem <- .rc_enrich_humangem_metadata_core(gem, reactions_tsv = reactions_tsv)
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
