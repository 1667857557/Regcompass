rc_layer2_penalty <- function(C_rel, Conf, epsilon = 1e-6, epsilon_C = 1e-3, epsilon_Conf = 1e-3,
                              penalty_cap = 20, support_reactions = character(), support_penalty = 0.05) {
  C_raw <- as.matrix(C_rel); Conf_raw <- as.matrix(Conf)
  Cc <- ifelse(is.na(C_raw), epsilon_C, pmax(C_raw, epsilon_C))
  Fc <- ifelse(is.na(Conf_raw), epsilon_Conf, pmax(Conf_raw, epsilon_Conf))
  E <- Cc * Fc
  P <- pmax(0, pmin(-log(E + epsilon), penalty_cap))
  P[!is.finite(P)] <- penalty_cap
  support_reaction_ids <- character()
  support_penalty_by_reaction <- numeric()
  if (length(support_reactions)) {
    if (is.null(names(support_reactions))) {
      support_reaction_ids <- intersect(as.character(support_reactions), rownames(P))
      support_penalty_by_reaction <- stats::setNames(rep(rc_layer2_support_penalty_for_type("support", support_penalty), length(support_reaction_ids)), support_reaction_ids)
    } else {
      support_reaction_ids <- intersect(names(support_reactions), rownames(P))
      support_penalty_by_reaction <- stats::setNames(as.numeric(support_reactions[support_reaction_ids]), support_reaction_ids)
    }
    if (length(support_penalty_by_reaction)) P[names(support_penalty_by_reaction), ] <- support_penalty_by_reaction
  }
  missing_evidence <- is.na(C_raw) | is.na(Conf_raw)
  list(penalty = P, components = list(C_rel = C_rel, reaction_confidence = Conf, evidence_product = E,
                                      missing_evidence = missing_evidence, support_reaction = rownames(P) %in% support_reaction_ids,
                                      support_penalty_by_reaction = support_penalty_by_reaction,
                                      penalty = P, epsilon_C = epsilon_C, epsilon_Conf = epsilon_Conf,
                                      penalty_cap = penalty_cap, support_penalty = support_penalty))
}
rc_compass_score_from_penalty <- function(P, feasible, epsilon = 1e-6) {
  score <- P
  for (i in seq_len(nrow(P))) {
    x <- P[i, ]; med <- stats::median(x[feasible[i, ]], na.rm = TRUE); sc <- stats::mad(x[feasible[i, ]], constant = 1.4826, na.rm = TRUE)
    if (!is.finite(sc) || sc <= 0) sc <- epsilon
    z <- (med - x) / (sc + epsilon); score[i, ] <- rc_sigmoid(z)
  }
  score[!feasible] <- NA_real_
  score[!is.finite(score) & feasible] <- NA_real_
  score
}

rc_layer2_unit_matrices <- function(layer1, unit, sample_col, celltype_col, condition_col) {
  C <- as.matrix(layer1$C_rel); Conf <- rc_layer2_confidence_matrix(layer1$reaction_confidence, C)
  common <- intersect(rownames(C), rownames(Conf)); C <- C[common,,drop=FALSE]; Conf <- Conf[common,,drop=FALSE]
  if (unit == "metacell") return(list(C_rel = C, reaction_confidence = Conf, unit_meta = layer1$unit_meta, summary = "metacell-level Layer1 matrices"))
  if (is.null(layer1$unit_meta)) stop("sample_celltype units require `layer1$unit_meta`.", call. = FALSE)
  pm <- layer1$unit_meta; required <- c("pool_id", sample_col, celltype_col); missing <- setdiff(required, colnames(pm)); if (length(missing)) stop("unit_meta missing: ", paste(missing, collapse=", "), call.=FALSE)
  pm <- pm[pm$pool_id %in% colnames(C), , drop=FALSE]
  group_cols <- c(sample_col, if (condition_col %in% colnames(pm)) condition_col else NULL, celltype_col)
  gid <- interaction(pm[, group_cols, drop=FALSE], sep="|", drop=TRUE)
  agg <- function(M) do.call(cbind, lapply(split(pm$pool_id, gid), function(cols) rowMedians_safe(M[, cols, drop=FALSE])))
  outC <- agg(C); outF <- agg(Conf)
  unit_meta <- unique(data.frame(unit_id = as.character(gid), pm[, group_cols, drop=FALSE], stringsAsFactors=FALSE))
  unit_meta <- unit_meta[match(colnames(outC), unit_meta$unit_id), , drop = FALSE]
  rownames(unit_meta) <- NULL; list(C_rel = outC, reaction_confidence = outF, unit_meta = unit_meta, summary = "Layer1 aggregated to sample × celltype medians")
}

rc_align_layer2_evidence <- function(M, rxns, fill = 1) {
  M <- as.matrix(M)
  out <- matrix(fill, nrow = length(rxns), ncol = ncol(M), dimnames = list(rxns, colnames(M)))
  common <- intersect(rxns, rownames(M))
  out[common, ] <- M[common, , drop = FALSE]
  out
}

rc_layer2_support_penalties <- function(gem, rxns, C_rel = NULL, Conf = NULL,
                                        support_classes = c("exchange", "demand", "sink", "support"),
                                        transport_penalty_mode = c("normal", "support", "reduced"),
                                        support_penalty = c(exchange = 0.05, demand = 0.1, sink = 0.1, support = 0.05),
                                        transport_reduced_penalty = 1) {
  transport_penalty_mode <- match.arg(transport_penalty_mode)
  meta <- gem$reaction_meta
  if (is.null(meta) || !is.data.frame(meta) || !"reaction_id" %in% colnames(meta)) return(numeric())
  type <- rc_layer2_reaction_type(meta)
  names(type) <- as.character(meta$reaction_id)
  type <- type[intersect(rxns, names(type))]
  out <- numeric()
  support_type <- type[type %in% support_classes]
  if (length(support_type)) {
    out <- c(out, stats::setNames(vapply(support_type, rc_layer2_support_penalty_for_type, numeric(1), support_penalty), names(support_type)))
  }
  transport <- type[type == "transport"]
  if (length(transport) && identical(transport_penalty_mode, "support")) {
    out <- c(out, stats::setNames(rep(rc_layer2_support_penalty_for_type("support", support_penalty), length(transport)), names(transport)))
  } else if (length(transport) && identical(transport_penalty_mode, "reduced")) {
    no_gpr <- !rc_layer2_has_gpr(meta[match(names(transport), as.character(meta$reaction_id)), , drop = FALSE], C_rel, Conf)
    if (any(no_gpr)) out <- c(out, stats::setNames(rep(transport_reduced_penalty, sum(no_gpr)), names(transport)[no_gpr]))
  }
  out[!duplicated(names(out))]
}

rc_layer2_reaction_type <- function(meta) {
  text_cols <- intersect(c("type", "reaction_type", "category", "subsystem", "name", "reaction_name"), colnames(meta))
  txt <- if (length(text_cols)) apply(meta[, text_cols, drop = FALSE], 1, paste, collapse = " ") else rep("", nrow(meta))
  txt <- tolower(txt)
  out <- rep("other", length(txt))
  out[out == "other" & grepl("exchange", txt)] <- "exchange"
  out[out == "other" & grepl("demand", txt)] <- "demand"
  out[out == "other" & grepl("sink", txt)] <- "sink"
  out[out == "other" & grepl("support|artificial", txt)] <- "support"
  out[out == "other" & grepl("transport", txt)] <- "transport"
  out
}

rc_layer2_support_penalty_for_type <- function(type, support_penalty) {
  if (length(support_penalty) == 1L || is.null(names(support_penalty))) return(as.numeric(support_penalty)[1])
  if (type %in% names(support_penalty)) return(as.numeric(support_penalty[[type]]))
  if ("support" %in% names(support_penalty)) return(as.numeric(support_penalty[["support"]]))
  as.numeric(support_penalty[[1]])
}

rc_layer2_has_gpr <- function(meta, C_rel = NULL, Conf = NULL) {
  gpr_cols <- intersect(c("gpr", "grRule", "gene_reaction_rule", "gene_reaction_rule_string", "genes"), colnames(meta))
  has_meta_gpr <- if (length(gpr_cols)) {
    apply(meta[, gpr_cols, drop = FALSE], 1, function(x) any(!is.na(x) & nzchar(trimws(as.character(x)))))
  } else {
    rep(FALSE, nrow(meta))
  }
  if (any(has_meta_gpr) || is.null(C_rel) || is.null(Conf)) return(has_meta_gpr)
  rid <- as.character(meta$reaction_id)
  C <- as.matrix(C_rel); F <- as.matrix(Conf)
  in_evidence <- rid %in% rownames(C) & rid %in% rownames(F)
  has_evidence <- rep(FALSE, length(rid))
  has_evidence[in_evidence] <- rowSums(is.finite(C[rid[in_evidence], , drop = FALSE]) | is.finite(F[rid[in_evidence], , drop = FALSE])) > 0
  has_evidence
}
rc_layer2_confidence_matrix <- function(confidence, C_rel) {
  C_rel <- as.matrix(C_rel)
  if (is.data.frame(confidence) && all(c("reaction_id", "pool_id", "reaction_confidence") %in% colnames(confidence))) {
    out <- matrix(NA_real_, nrow = nrow(C_rel), ncol = ncol(C_rel), dimnames = dimnames(C_rel))
    rid <- as.character(confidence$reaction_id)
    pid <- as.character(confidence$pool_id)
    ok <- rid %in% rownames(out) & pid %in% colnames(out)
    if (any(ok)) {
      idx <- cbind(match(rid[ok], rownames(out)), match(pid[ok], colnames(out)))
      out[idx] <- as.numeric(confidence$reaction_confidence[ok])
    }
    return(out)
  }
  M <- as.matrix(confidence)
  if (is.null(rownames(M)) || is.null(colnames(M))) stop("`reaction_confidence` must be a matrix-like object with dimnames or a long data frame.", call. = FALSE)
  M
}
