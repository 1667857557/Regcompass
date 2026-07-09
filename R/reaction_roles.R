
.rc_meta_col <- function(meta, names, n) {
  nm <- intersect(names, colnames(meta))[1]
  if (is.na(nm)) return(rep("", n))
  x <- as.character(meta[[nm]])
  x[is.na(x)] <- ""
  tolower(x)
}

.rc_infer_exchange_like <- function(gem, gv, meta, role) {
  n <- length(gv$reactions)
  rid <- tolower(as.character(gv$reactions))
  rn <- .rc_meta_col(meta, c("reaction_name", "name", "description"), n)
  rs <- .rc_meta_col(meta, c("subsystem", "subSystems", "subSystem", "sub_system", "metabolic_module"), n)
  eq <- .rc_meta_col(meta, c("equation", "reaction_formula", "formula"), n)

  id_exchange <- grepl("^(ex_|ex|exchange)", rid)
  name_exchange <- grepl("exchange|boundary exchange|extracellular exchange|uptake|secretion", rn)
  subsystem_exchange <- grepl("exchange|extracellular", rs)
  eq_exchange <- grepl("\\[e\\]|\\[extracellular\\]|extracellular", eq)

  comp <- NULL
  if (!is.null(gem$metabolite_meta) && "compartment" %in% colnames(gem$metabolite_meta)) {
    comp <- as.character(gem$metabolite_meta$compartment[
      match(rownames(gv$S), as.character(gem$metabolite_meta$metabolite_id))
    ])
  }

  single_extracellular <- rep(FALSE, n)
  if (!is.null(comp)) {
    single_extracellular <- vapply(seq_len(ncol(gv$S)), function(j) {
      idx <- which(gv$S[, j] != 0)
      rxn_comp <- unique(stats::na.omit(comp[idx]))
      length(idx) == 1L && length(rxn_comp) == 1L && identical(rxn_comp[[1]], "e")
    }, logical(1))
  }

  role == "unknown" & (
    id_exchange |
      name_exchange |
      subsystem_exchange |
      eq_exchange |
      single_extracellular
  )
}

#' Annotate GEM reactions with curated or inferred RegCompassR roles
#' @export
rc_annotate_reaction_roles <- function(gem, reaction_role_table = NULL, infer_from_id = TRUE,
                                       infer_from_stoichiometry = TRUE, infer_from_compartment = TRUE,
                                       overwrite_existing = FALSE) {
  gv <- rc_validate_gem(gem)
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta)) meta <- data.frame(reaction_id = gv$reactions, stringsAsFactors = FALSE)
  if (!"reaction_id" %in% colnames(meta)) meta$reaction_id <- gv$reactions
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  meta$reaction_id <- gv$reactions
  old_role <- if ("role" %in% colnames(meta)) as.character(meta$role) else rep(NA_character_, nrow(meta))
  role <- ifelse(!is.na(old_role) & nzchar(old_role) & !overwrite_existing, old_role, "unknown")
  source <- ifelse(role != "unknown", "metadata", "unknown")
  conf <- ifelse(role != "unknown", "medium", "low")
  if (infer_from_id || infer_from_stoichiometry || infer_from_compartment) {
    idx_ex <- .rc_infer_exchange_like(gem, gv, meta, role)
    role[idx_ex] <- "exchange"; source[idx_ex] <- "metadata_stoichiometry"; conf[idx_ex] <- "medium"
  }
  if (infer_from_stoichiometry) {
    nnz <- Matrix::colSums(abs(gv$S) > 0)
    idx <- role == "unknown" & nnz == 1
    role[idx] <- "boundary_like"; source[idx] <- "stoichiometry"; conf[idx] <- "low"
  }
  if (infer_from_compartment && !is.null(gem$metabolite_meta) && "compartment" %in% colnames(gem$metabolite_meta)) {
    comp <- as.character(gem$metabolite_meta$compartment[match(rownames(gv$S), as.character(gem$metabolite_meta$metabolite_id))])
    comps_by_rxn <- lapply(seq_len(ncol(gv$S)), function(j) unique(stats::na.omit(comp[which(gv$S[, j] != 0)])))
    idx <- role == "unknown" & vapply(comps_by_rxn, length, integer(1)) > 1L
    role[idx] <- "transport"; source[idx] <- "stoichiometry"; conf[idx] <- "medium"
  }
  if (infer_from_id) {
    rid <- tolower(gv$reactions)
    pats <- list(exchange = "^(ex_|ex)", demand = "^(dm_|demand)", sink = "^(sink_|sk_)", biomass = "biomass", maintenance = "maintenance|atpm", transport = "(_tx$|transport|t$)")
    for (nm in names(pats)) {
      idx <- role == "unknown" & grepl(pats[[nm]], rid)
      role[idx] <- nm; source[idx] <- "id_pattern"; conf[idx] <- "low"
    }
  }
  if (!is.null(reaction_role_table)) {
    if (!all(c("reaction_id", "role") %in% colnames(reaction_role_table))) stop("`reaction_role_table` must contain `reaction_id` and `role`.", call. = FALSE)
    m <- match(gv$reactions, as.character(reaction_role_table$reaction_id))
    idx <- !is.na(m)
    role[idx] <- as.character(reaction_role_table$role[m[idx]])
    source[idx] <- "curated"; conf[idx] <- "high"
  }
  meta$role <- role; meta$role_source <- source; meta$role_confidence <- conf
  gem$reaction_meta <- meta; gem$reaction_roles <- meta[, c("reaction_id", "role", "role_source", "role_confidence"), drop = FALSE]
  gem
}
