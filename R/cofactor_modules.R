#' Extract dynamic cofactor-regeneration module candidates from GEM annotations
#' @export
rc_extract_cofactor_modules <- function(gem,
                                        cofactors = c("nadph", "nadp", "nadh", "nad", "atp", "adp", "coa"),
                                        subsystem_patterns = c("pentose phosphate", "folate", "one carbon", "TCA", "oxidative phosphorylation", "malic enzyme", "isocitrate", "glutathione", "fatty acid"),
                                        require_gpr = TRUE) {
  gv <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  if (is.null(meta)) meta <- data.frame(reaction_id = gv$reactions, stringsAsFactors = FALSE)
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  subsystem <- if ("subsystem" %in% colnames(meta)) as.character(meta$subsystem) else if ("subSystems" %in% colnames(meta)) as.character(meta$subSystems) else rep("", ncol(gv$S))
  role <- if ("role" %in% colnames(meta)) as.character(meta$role) else rep(NA_character_, ncol(gv$S))
  gpr <- if ("gpr" %in% colnames(meta)) as.character(meta$gpr) else if ("grRules" %in% colnames(meta)) as.character(meta$grRules) else rep("", ncol(gv$S))
  has_gpr <- nzchar(gpr) | gv$reactions %in% unique(as.character((gem$gpr_table %||% data.frame(reaction_id=character()))$reaction_id))
  met_ids <- tolower(rownames(gv$S))
  touches <- vapply(seq_len(ncol(gv$S)), function(j) any(vapply(cofactors, function(cf) any(grepl(cf, met_ids[gv$S[, j] != 0], fixed = TRUE)), logical(1))), logical(1))
  pat <- paste(subsystem_patterns, collapse = "|")
  sub_hit <- grepl(pat, subsystem, ignore.case = TRUE)
  role_hit <- grepl("cofactor|redox|energy|transport", role, ignore.case = TRUE)
  keep <- touches & (sub_hit | role_hit | has_gpr)
  if (isTRUE(require_gpr)) keep <- keep & (has_gpr | role_hit)
  rxn_cof <- lapply(which(keep), function(j) cofactors[vapply(cofactors, function(cf) any(grepl(cf, met_ids[gv$S[, j] != 0], fixed = TRUE)), logical(1))])
  pairs <- vapply(rxn_cof, paste, collapse = "/", FUN.VALUE = character(1))
  data.frame(cofactor_pair = pairs,
             candidate_reactions = gv$reactions[keep],
             selected_reactions = gv$reactions[keep],
             selection_reason = ifelse(sub_hit[keep], "subsystem_pattern", ifelse(role_hit[keep], "reaction_role", "gpr_cofactor_stoichiometry")),
             has_gpr = has_gpr[keep],
             median_C_rel = NA_real_,
             median_confidence = NA_real_,
             cofactor_limitation_flag = NA,
             stringsAsFactors = FALSE)
}
