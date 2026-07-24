# Biological reaction annotations and evidence joins.

.rc_ra_nonempty <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(trimws(x))
}

.rc_ra_first_col <- function(x, candidates) {
  hit <- candidates[candidates %in% colnames(x)]
  if (length(hit)) hit[[1L]] else NULL
}

.rc_ra_collapse <- function(x, sep = ";") {
  x <- trimws(as.character(x))
  x <- unique(x[.rc_ra_nonempty(x)])
  if (length(x)) paste(x, collapse = sep) else NA_character_
}

.rc_ra_split <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x))) return(character())
  out <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  unique(out[nzchar(out)])
}

.rc_ra_format_coefficient <- function(x, tolerance = 1e-12) {
  x <- abs(as.numeric(x))
  if (!is.finite(x) || abs(x - 1) <= tolerance) return("")
  paste0(format(signif(x, 8), trim = TRUE, scientific = FALSE), " ")
}

.rc_ra_metabolite_labels <- function(gem) {
  ids <- rownames(gem$S)
  meta <- gem$metabolite_meta
  if (is.null(meta) || !is.data.frame(meta)) {
    return(stats::setNames(ids, ids))
  }
  id_col <- .rc_ra_first_col(meta, c("metabolite_id", "mets", "id"))
  if (is.null(id_col)) return(stats::setNames(ids, ids))
  idx <- match(ids, as.character(meta[[id_col]]))
  name_col <- .rc_ra_first_col(meta, c("name", "metNames", "metabolite_name"))
  comp_col <- .rc_ra_first_col(meta, c("compartment", "compartments", "metComps"))
  labels <- ids
  if (!is.null(name_col)) {
    names_in <- as.character(meta[[name_col]][idx])
    use <- .rc_ra_nonempty(names_in)
    labels[use] <- names_in[use]
  }
  if (!is.null(comp_col)) {
    compartments <- as.character(meta[[comp_col]][idx])
    use <- .rc_ra_nonempty(compartments)
    labels[use] <- paste0(labels[use], " [", compartments[use], "]")
  }
  stats::setNames(labels, ids)
}

.rc_ra_formula_parts <- function(stoichiometry, labels) {
  negative <- which(is.finite(stoichiometry) & stoichiometry < 0)
  positive <- which(is.finite(stoichiometry) & stoichiometry > 0)
  make_side <- function(index) {
    if (!length(index)) return("empty")
    paste(vapply(index, function(i) {
      paste0(
        .rc_ra_format_coefficient(stoichiometry[[i]]),
        labels[[i]]
      )
    }, character(1)), collapse = " + ")
  }
  list(left = make_side(negative), right = make_side(positive))
}

.rc_ra_gpr_from_table <- function(gem, layer1 = NULL) {
  table <- gem$gpr_table
  if (!is.null(table) && is.data.frame(table) &&
      all(c("reaction_id", "and_group_id", "gene") %in% colnames(table))) {
    table <- table[.rc_ra_nonempty(table$reaction_id) &
                     .rc_ra_nonempty(table$gene), , drop = FALSE]
    if (nrow(table)) {
      table$reaction_id <- as.character(table$reaction_id)
      table$gene <- toupper(trimws(as.character(table$gene)))
      return(split(table, table$reaction_id))
    }
  }
  parsed <- if (is.list(layer1)) layer1$parsed_gpr else NULL
  if (!is.list(parsed) || !length(parsed)) return(list())
  lapply(names(parsed), function(reaction_id) {
    groups <- parsed[[reaction_id]]
    rows <- lapply(seq_along(groups), function(i) {
      genes <- toupper(trimws(as.character(groups[[i]])))
      genes <- unique(genes[.rc_ra_nonempty(genes)])
      if (!length(genes)) return(NULL)
      data.frame(
        reaction_id = reaction_id,
        and_group_id = i,
        gene = genes,
        stringsAsFactors = FALSE
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) do.call(rbind, rows) else NULL
  }) |>
    stats::setNames(names(parsed))
}

.rc_ra_gpr_strings <- function(gpr_rows) {
  if (is.null(gpr_rows) || !is.data.frame(gpr_rows) || !nrow(gpr_rows)) {
    return(list(genes = NA_character_, gpr_rule = NA_character_, n_genes = 0L))
  }
  genes <- sort(unique(toupper(trimws(as.character(gpr_rows$gene)))))
  genes <- genes[.rc_ra_nonempty(genes)]
  groups <- split(gpr_rows, as.character(gpr_rows$and_group_id))
  group_rules <- vapply(groups, function(one) {
    values <- unique(toupper(trimws(as.character(one$gene))))
    values <- values[.rc_ra_nonempty(values)]
    if (!length(values)) return(NA_character_)
    rule <- paste(values, collapse = " and ")
    if (length(values) > 1L) paste0("(", rule, ")") else rule
  }, character(1))
  group_rules <- group_rules[.rc_ra_nonempty(group_rules)]
  list(
    genes = if (length(genes)) paste(genes, collapse = ";") else NA_character_,
    gpr_rule = if (length(group_rules)) paste(group_rules, collapse = " or ") else NA_character_,
    n_genes = length(genes)
  )
}

.rc_ra_reaction_catalog <- function(gem, layer1 = NULL, reaction_ids = NULL) {
  rc_validate_gem(gem)
  S <- gem$S
  all_ids <- colnames(S)
  if (is.null(reaction_ids)) reaction_ids <- all_ids
  reaction_ids <- unique(as.character(reaction_ids))
  unknown <- setdiff(reaction_ids, all_ids)
  if (length(unknown)) {
    stop(
      "Reaction IDs are absent from the GEM: ",
      paste(utils::head(unknown, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  reaction_ids <- all_ids[all_ids %in% reaction_ids]
  labels <- .rc_ra_metabolite_labels(gem)
  reaction_meta <- gem$reaction_meta
  if (is.null(reaction_meta) || !is.data.frame(reaction_meta)) {
    reaction_meta <- data.frame(reaction_id = all_ids, stringsAsFactors = FALSE)
  }
  id_col <- .rc_ra_first_col(reaction_meta, c("reaction_id", "rxns", "id"))
  if (is.null(id_col)) {
    reaction_meta$reaction_id <- all_ids
    id_col <- "reaction_id"
  }
  meta_idx <- match(reaction_ids, as.character(reaction_meta[[id_col]]))
  name_col <- .rc_ra_first_col(reaction_meta, c("name", "rxnNames", "reaction_name"))
  subsystem_col <- .rc_ra_first_col(
    reaction_meta, c("subsystem", "subSystems", "metabolic_module")
  )
  role_col <- .rc_ra_first_col(reaction_meta, c("role", "reaction_role"))
  crossref_candidates <- list(
    kegg_reaction_id = c("kegg_reaction_id", "rxnKEGGID", "kegg_id"),
    reactome_reaction_id = c("reactome_reaction_id", "rxnREACTOMEID", "reactome_id"),
    rhea_reaction_id = c("rhea_reaction_id", "rxnRheaID", "rhea_id"),
    rhea_master_id = c("rhea_master_id", "rxnRheaMasterID", "master_rhea_id")
  )
  gpr_split <- .rc_ra_gpr_from_table(gem, layer1)
  rows <- lapply(reaction_ids, function(reaction_id) {
    j <- match(reaction_id, all_ids)
    stoichiometry <- as.numeric(S[, j, drop = TRUE])
    names(stoichiometry) <- rownames(S)
    nonzero <- which(is.finite(stoichiometry) & abs(stoichiometry) > 1e-12)
    parts <- .rc_ra_formula_parts(
      stoichiometry[nonzero],
      unname(labels[names(stoichiometry)[nonzero]])
    )
    lb <- as.numeric(gem$lb[[reaction_id]])
    ub <- as.numeric(gem$ub[[reaction_id]])
    arrow <- if (is.finite(lb) && is.finite(ub) && lb < 0 && ub > 0) {
      " <=> "
    } else if (is.finite(ub) && ub > 0) {
      " -> "
    } else if (is.finite(lb) && lb < 0) {
      " <- "
    } else {
      " -/-> "
    }
    idx <- meta_idx[match(reaction_id, reaction_ids)]
    reaction_name <- if (!is.null(name_col) && is.finite(idx)) {
      as.character(reaction_meta[[name_col]][idx])
    } else {
      NA_character_
    }
    if (!.rc_ra_nonempty(reaction_name)) reaction_name <- reaction_id
    subsystem <- if (!is.null(subsystem_col) && is.finite(idx)) {
      as.character(reaction_meta[[subsystem_col]][idx])
    } else {
      NA_character_
    }
    role <- if (!is.null(role_col) && is.finite(idx)) {
      as.character(reaction_meta[[role_col]][idx])
    } else {
      "internal"
    }
    if (!.rc_ra_nonempty(role)) role <- "internal"
    gpr <- .rc_ra_gpr_strings(gpr_split[[reaction_id]])
    row <- data.frame(
      reaction_id = reaction_id,
      reaction_name = reaction_name,
      subsystem = subsystem,
      reaction_role = role,
      lower_bound = lb,
      upper_bound = ub,
      reversible = is.finite(lb) && is.finite(ub) && lb < 0 && ub > 0,
      model_formula = paste0(parts$left, arrow, parts$right),
      forward_substrates = parts$left,
      forward_products = parts$right,
      forward_formula = paste0(parts$left, " -> ", parts$right),
      reverse_substrates = parts$right,
      reverse_products = parts$left,
      reverse_formula = paste0(parts$right, " -> ", parts$left),
      genes = gpr$genes,
      gpr_rule = gpr$gpr_rule,
      n_gpr_genes = gpr$n_genes,
      has_gpr = gpr$n_genes > 0L,
      stringsAsFactors = FALSE
    )
    for (target in names(crossref_candidates)) {
      column <- .rc_ra_first_col(reaction_meta, crossref_candidates[[target]])
      row[[target]] <- if (!is.null(column) && is.finite(idx)) {
        as.character(reaction_meta[[column]][idx])
      } else {
        NA_character_
      }
    }
    row
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.rc_ra_infer_group_columns <- function(meta, condition_col, celltype_col) {
  if (!is.data.frame(meta)) {
    stop("Layer 1 `unit_meta` must be a data frame.", call. = FALSE)
  }
  if (is.null(condition_col)) {
    condition_col <- .rc_ra_first_col(
      meta, c("condition", "dataset", "Group", "group", "treatment")
    )
  }
  if (is.null(celltype_col)) {
    celltype_col <- .rc_ra_first_col(
      meta, c("cell_type", "celltype", "epithelial_or_stem", "CellType")
    )
  }
  if (is.null(condition_col) || is.null(celltype_col) ||
      !condition_col %in% colnames(meta) || !celltype_col %in% colnames(meta)) {
    stop(
      "Could not identify condition and cell-type columns in Layer 1 metadata.",
      call. = FALSE
    )
  }
  list(condition_col = condition_col, celltype_col = celltype_col)
}

.rc_ra_annotation_from_object <- function(x) {
  if (!is.list(x)) return(NULL)
  catalog <- x$reaction_catalog
  evidence <- x$reaction_evidence
  if (is.data.frame(catalog)) {
    return(list(
      reactions = catalog,
      evidence = if (is.data.frame(evidence)) evidence else data.frame()
    ))
  }
  NULL
}

.rc_ra_catalog_join <- function(data, catalog) {
  if (!is.data.frame(data) || !nrow(data) ||
      !is.data.frame(catalog) || !nrow(catalog)) return(data)
  idx <- match(as.character(data$reaction_id), as.character(catalog$reaction_id))
  columns <- setdiff(colnames(catalog), "reaction_id")
  for (column in columns) data[[column]] <- catalog[[column]][idx]
  reverse <- as.character(data$target_direction) == "reverse"
  data$tested_substrates <- ifelse(
    reverse, data$reverse_substrates, data$forward_substrates
  )
  data$tested_products <- ifelse(
    reverse, data$reverse_products, data$forward_products
  )
  data$tested_formula <- ifelse(
    reverse, data$reverse_formula, data$forward_formula
  )
  data
}

.rc_ra_evidence_index <- function(evidence) {
  if (!is.data.frame(evidence) || !nrow(evidence)) return(character())
  paste(
    evidence$reaction_id,
    evidence$condition,
    evidence$cell_type,
    sep = "\001"
  )
}

.rc_ra_group_evidence_join <- function(data, evidence, condition_col = "condition") {
  if (!is.data.frame(data) || !nrow(data) ||
      !is.data.frame(evidence) || !nrow(evidence)) return(data)
  key <- paste(
    data$reaction_id,
    as.character(data[[condition_col]]),
    data$cell_type,
    sep = "\001"
  )
  idx <- match(key, .rc_ra_evidence_index(evidence))
  fields <- setdiff(
    colnames(evidence), c("reaction_id", "condition", "cell_type")
  )
  for (field in fields) data[[field]] <- evidence[[field]][idx]
  data
}

.rc_ra_enrich_statistics <- function(answer, annotation) {
  if (is.null(annotation)) return(answer)
  answer$pairwise <- .rc_ra_catalog_join(
    answer$pairwise, annotation$reactions
  )
  answer$pairwise <- .rc_ra_pairwise_evidence(
    answer$pairwise, annotation$evidence
  )
  answer$omnibus <- .rc_ra_catalog_join(
    answer$omnibus, annotation$reactions
  )
  answer$omnibus <- .rc_ra_omnibus_evidence(
    answer$omnibus,
    annotation$evidence,
    answer$params$conditions
  )
  reaction_ids <- unique(c(
    as.character(answer$pairwise$reaction_id),
    as.character(answer$omnibus$reaction_id)
  ))
  answer$reaction_catalog <- annotation$reactions[
    annotation$reactions$reaction_id %in% reaction_ids,
    , drop = FALSE
  ]
  answer$reaction_evidence <- annotation$evidence[
    annotation$evidence$reaction_id %in% reaction_ids,
    , drop = FALSE
  ]
  answer$params$reaction_annotations <- TRUE
  answer
}

.rc_ra_scored_reaction_ids <- function(layer2) {
  if (!is.list(layer2) || is.null(layer2$penalty) ||
      is.null(rownames(layer2$penalty))) return(NULL)
  unique(as.character(
    rc_parse_microcompass_row_id(rownames(layer2$penalty))$reaction_id
  ))
}

.rc_ra_attach_to_result <- function(
    result, gem, condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  if (!is.list(result) || is.null(result$layer1) || is.null(result$microcompass)) {
    stop("`result` must contain Layer 1 and microCOMPASS outputs.", call. = FALSE)
  }
  reaction_ids <- .rc_ra_scored_reaction_ids(result$microcompass)
  annotation <- rc_build_reaction_annotations(
    gem = gem,
    layer1 = result$layer1,
    reaction_ids = reaction_ids,
    condition_col = condition_col,
    celltype_col = celltype_col,
    evidence_tolerance = evidence_tolerance
  )
  result$reaction_catalog <- annotation$reactions
  result$reaction_evidence <- annotation$evidence
  result$reaction_ranking <- .rc_ra_catalog_join(
    result$reaction_ranking, annotation$reactions
  )
  result$reaction_ranking <- .rc_ra_group_evidence_join(
    result$reaction_ranking, annotation$evidence, "condition"
  )
  result$condition_summary <- .rc_ra_catalog_join(
    result$condition_summary, annotation$reactions
  )
  result$condition_summary <- .rc_ra_group_evidence_join(
    result$condition_summary, annotation$evidence, "condition"
  )
  result$condition_contrast <- .rc_ra_catalog_join(
    result$condition_contrast, annotation$reactions
  )
  result$condition_contrast <- .rc_ra_pairwise_evidence(
    result$condition_contrast, annotation$evidence
  )
  result$params$reaction_annotation <- list(
    available = TRUE,
    n_reactions = nrow(annotation$reactions),
    n_group_evidence_rows = nrow(annotation$evidence),
    evidence_tolerance = evidence_tolerance,
    evidence_classes = c(
      "RNA+ATAC", "RNA-only", "GPR/no-observed-RNA", "structural/no-GPR"
    )
  )
  result
}

#' Attach reaction annotations to an existing RegCompass result
#'
#' Adds formal reaction names, formulas, GPR genes, and condition-by-cell-type
#' RNA versus RNA+ATAC evidence provenance to a result produced by an earlier
#' RegCompass version.
#'
#' @param result Complete RegCompass result.
#' @param gem The exact GEM used for the analysis.
#' @param condition_col,celltype_col Optional Layer 1 metadata columns.
#' @param evidence_tolerance Numerical evidence tolerance.
#' @return The annotated RegCompass result.
#' @export
rc_attach_reaction_annotations <- function(
    result, gem, condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  .rc_ra_attach_to_result(
    result = result,
    gem = gem,
    condition_col = condition_col,
    celltype_col = celltype_col,
    evidence_tolerance = evidence_tolerance
  )
}

#' Select scored reactions by metabolic genes
#'
#' Selects reactions whose GPR contains any or all requested genes. The returned
#' reaction table includes formal names, formulas, GPR rules, and the matched
#' genes. Evidence rows can be restricted to selected conditions, cell types, or
#' evidence classes.
#'
#' @param x An annotated RegCompass result.
#' @param genes Metabolic gene symbols.
#' @param match Require `"any"` or `"all"` requested genes in each reaction.
#' @param conditions,cell_types Optional evidence-group filters.
#' @param evidence_class Optional evidence classes such as `"RNA+ATAC"` or
#'   `"RNA-only"`.
#' @return A `regcompass_gene_reaction_selection` list.
#' @export
rc_select_gene_reactions <- function(
    x, genes, match = c("any", "all"),
    conditions = NULL, cell_types = NULL, evidence_class = NULL) {
  match <- match.arg(match)
  annotation <- .rc_ra_annotation_from_object(x)
  if (is.null(annotation)) {
    stop(
      "Reaction annotations are unavailable. Re-run Stage 6 or call ",
      "`rc_attach_reaction_annotations(result, gem)` first.",
      call. = FALSE
    )
  }
  genes <- unique(toupper(trimws(as.character(genes))))
  genes <- genes[.rc_ra_nonempty(genes)]
  if (!length(genes)) stop("`genes` must contain at least one symbol.", call. = FALSE)
  catalog <- annotation$reactions
  matched <- lapply(catalog$genes, function(value) {
    intersect(genes, toupper(.rc_ra_split(value)))
  })
  keep <- if (identical(match, "all")) {
    vapply(matched, length, integer(1)) == length(genes)
  } else {
    vapply(matched, length, integer(1)) > 0L
  }
  catalog <- catalog[keep, , drop = FALSE]
  matched <- matched[keep]
  catalog$matched_genes <- vapply(matched, .rc_ra_collapse, character(1))
  catalog$n_matched_genes <- vapply(matched, length, integer(1))
  evidence <- annotation$evidence
  if (is.data.frame(evidence) && nrow(evidence)) {
    evidence <- evidence[evidence$reaction_id %in% catalog$reaction_id, , drop = FALSE]
    if (!is.null(conditions)) {
      evidence <- evidence[evidence$condition %in% as.character(conditions), , drop = FALSE]
    }
    if (!is.null(cell_types)) {
      evidence <- evidence[evidence$cell_type %in% as.character(cell_types), , drop = FALSE]
    }
    if (!is.null(evidence_class)) {
      evidence <- evidence[
        evidence$evidence_class %in% as.character(evidence_class), , drop = FALSE
      ]
      catalog <- catalog[
        catalog$reaction_id %in% evidence$reaction_id, , drop = FALSE
      ]
    }
  }
  targets <- data.frame()
  if (is.list(x) && is.list(x$microcompass) &&
      !is.null(x$microcompass$penalty)) {
    targets <- rc_parse_microcompass_row_id(rownames(x$microcompass$penalty))
    targets$row_id <- rownames(x$microcompass$penalty)
    targets <- targets[targets$reaction_id %in% catalog$reaction_id, , drop = FALSE]
    targets <- .rc_ra_catalog_join(targets, catalog)
  }
  answer <- list(
    reactions = catalog,
    evidence = evidence,
    targets = targets,
    reaction_ids = unique(as.character(catalog$reaction_id)),
    genes = genes,
    match = match
  )
  class(answer) <- c("regcompass_gene_reaction_selection", "list")
  answer
}

.rc_ra_plot_caption <- function(annotation_row, evidence_text = NULL) {
  if (!is.data.frame(annotation_row) || !nrow(annotation_row)) return(NULL)
  parts <- c(
    if (.rc_ra_nonempty(annotation_row$tested_formula[[1L]])) {
      annotation_row$tested_formula[[1L]]
    },
    if (.rc_ra_nonempty(annotation_row$genes[[1L]])) {
      paste0("Genes: ", annotation_row$genes[[1L]])
    },
    if (.rc_ra_nonempty(evidence_text)) {
      paste0("Evidence: ", evidence_text)
    }
  )
  if (length(parts)) paste(parts, collapse = "\n") else NULL
}



.rc_ra_annotate_condition_plot <- function(
    plot, statistics, reaction_id, cell_type,
    target_direction = NULL, medium_scenario = NULL, title = NULL) {
  annotation_row <- data.frame()
  if (!is.list(statistics) || !is.data.frame(statistics$pairwise)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  required <- c(
    "reaction_id", "cell_type", "target_direction", "medium_scenario",
    "reaction_name", "tested_formula", "genes", "evidence_comparison"
  )
  if (!all(required %in% colnames(statistics$pairwise))) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }

  keep <- as.character(statistics$pairwise$reaction_id) == reaction_id &
    as.character(statistics$pairwise$cell_type) == cell_type
  if (!is.null(target_direction)) {
    keep <- keep &
      as.character(statistics$pairwise$target_direction) == target_direction
  }
  if (!is.null(medium_scenario)) {
    keep <- keep &
      as.character(statistics$pairwise$medium_scenario) == medium_scenario
  }
  annotation_row <- statistics$pairwise[keep, , drop = FALSE]
  if (!nrow(annotation_row)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  annotation_row <- annotation_row[1L, , drop = FALSE]
  reaction_name <- as.character(annotation_row$reaction_name[[1L]])
  evidence_text <- as.character(annotation_row$evidence_comparison[[1L]])
  if (is.null(title) && length(reaction_name) == 1L &&
      .rc_ra_nonempty(reaction_name)) {
    plot <- plot + ggplot2::labs(
      title = paste0(
        reaction_name, " (", reaction_id, ", ",
        annotation_row$target_direction[[1L]], ") in ", cell_type
      )
    )
  }
  caption <- .rc_ra_plot_caption(annotation_row, evidence_text)
  if (!is.null(caption) && length(caption) == 1L && !is.na(caption)) {
    plot <- plot + ggplot2::labs(caption = caption) + ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  }
  attr(plot, "reaction_annotation") <- annotation_row
  plot
}
