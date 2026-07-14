.rc_humangem_met_compartment_from_id <- function(metabolite_id) {
  x <- as.character(metabolite_id)
  known <- c("c", "e", "g", "l", "m", "n", "p", "r", "x")
  bracket <- sub("^.*\\[([A-Za-z])\\]$", "\\1", x)
  last <- substring(x, nchar(x), nchar(x))
  out <- ifelse(bracket %in% known, bracket, ifelse(last %in% known, last, NA_character_))
  out
}

#' Enrich Human-GEM metadata with module, equation, reversibility, and compartment annotations
#' @export
rc_enrich_humangem_metadata <- function(gem, reactions_tsv = NULL, metabolites_tsv = NULL) {
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
