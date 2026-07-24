# Biological reaction annotations and gene-centered condition analysis.
# This file is collated after the canonical statistics, plotting, and result
# assembly implementations so the public entry points can enrich their outputs
# without changing the underlying LP or statistical calculations.

.rc_condition_statistics_core <- rc_test_condition_reactions
.rc_condition_plot_core <- rc_plot_condition_reaction
.rc_step_results_core <- rc_regcompass_step_results

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

.rc_ra_group_evidence <- function(
    catalog, layer1, condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  empty <- data.frame(
    reaction_id = character(), condition = character(), cell_type = character(),
    evidence_class = character(), n_units = integer(),
    rna_supported_genes = character(), n_rna_supported_genes = integer(),
    atac_modifier_genes = character(), n_atac_modifier_genes = integer(),
    multiome_contributing_genes = character(),
    n_multiome_contributing_genes = integer(),
    has_rna_evidence = logical(),
    has_atac_regulatory_evidence = logical(),
    has_active_multiome_contribution = logical(),
    stringsAsFactors = FALSE
  )
  if (is.null(layer1) || !is.list(layer1)) return(empty)
  required <- c(
    "gene_support_rna", "gene_regulatory_modifier",
    "gene_support_multiome", "unit_meta"
  )
  if (!all(required %in% names(layer1))) return(empty)
  rna <- as.matrix(layer1$gene_support_rna)
  modifier <- as.matrix(layer1$gene_regulatory_modifier)
  multiome <- as.matrix(layer1$gene_support_multiome)
  if (!identical(dimnames(rna), dimnames(modifier)) ||
      !identical(dimnames(rna), dimnames(multiome))) {
    stop("Layer 1 RNA, ATAC modifier, and multiome matrices must align.", call. = FALSE)
  }
  meta <- layer1$unit_meta
  columns <- .rc_ra_infer_group_columns(meta, condition_col, celltype_col)
  condition_col <- columns$condition_col
  celltype_col <- columns$celltype_col
  unit_col <- .rc_ra_first_col(meta, c("unit_id", "pool_id", "metacell_id"))
  if (is.null(unit_col)) stop("Layer 1 metadata lack unit identifiers.", call. = FALSE)
  unit_ids <- as.character(meta[[unit_col]])
  if (!setequal(colnames(rna), unit_ids)) {
    stop("Layer 1 evidence matrices and metadata contain different units.", call. = FALSE)
  }
  meta <- meta[match(colnames(rna), unit_ids), , drop = FALSE]
  rownames(rna) <- tolower(rownames(rna))
  rownames(modifier) <- tolower(rownames(modifier))
  rownames(multiome) <- tolower(rownames(multiome))
  groups <- unique(data.frame(
    condition = as.character(meta[[condition_col]]),
    cell_type = as.character(meta[[celltype_col]]),
    stringsAsFactors = FALSE
  ))
  groups <- groups[.rc_ra_nonempty(groups$condition) &
                     .rc_ra_nonempty(groups$cell_type), , drop = FALSE]
  gene_lists <- lapply(catalog$genes, function(x) {
    tolower(.rc_ra_split(x))
  })
  names(gene_lists) <- catalog$reaction_id
  output <- vector("list", nrow(groups))
  for (group_index in seq_len(nrow(groups))) {
    condition <- groups$condition[[group_index]]
    cell_type <- groups$cell_type[[group_index]]
    keep <- as.character(meta[[condition_col]]) == condition &
      as.character(meta[[celltype_col]]) == cell_type
    rna_group <- rna[, keep, drop = FALSE]
    modifier_group <- modifier[, keep, drop = FALSE]
    multiome_group <- multiome[, keep, drop = FALSE]
    rna_flag <- rowSums(
      is.finite(rna_group) & rna_group > evidence_tolerance,
      na.rm = TRUE
    ) > 0
    modifier_flag <- rowSums(
      is.finite(modifier_group) & abs(modifier_group) > evidence_tolerance,
      na.rm = TRUE
    ) > 0
    contribution_flag <- rowSums(
      is.finite(multiome_group) & is.finite(rna_group) &
        abs(multiome_group - rna_group) > evidence_tolerance,
      na.rm = TRUE
    ) > 0
    names(rna_flag) <- rownames(rna)
    names(modifier_flag) <- rownames(rna)
    names(contribution_flag) <- rownames(rna)
    rows <- lapply(catalog$reaction_id, function(reaction_id) {
      genes <- gene_lists[[reaction_id]]
      genes <- intersect(genes, rownames(rna))
      rna_genes <- genes[unname(rna_flag[genes])]
      modifier_genes <- genes[unname(modifier_flag[genes])]
      contribution_genes <- genes[unname(contribution_flag[genes])]
      evidence_class <- if (!length(.rc_ra_split(
        catalog$genes[match(reaction_id, catalog$reaction_id)]
      ))) {
        "structural/no-GPR"
      } else if (!length(rna_genes)) {
        "GPR/no-observed-RNA"
      } else if (length(contribution_genes)) {
        "RNA+ATAC"
      } else {
        "RNA-only"
      }
      data.frame(
        reaction_id = reaction_id,
        condition = condition,
        cell_type = cell_type,
        evidence_class = evidence_class,
        n_units = sum(keep),
        rna_supported_genes = .rc_ra_collapse(toupper(rna_genes)),
        n_rna_supported_genes = length(rna_genes),
        atac_modifier_genes = .rc_ra_collapse(toupper(modifier_genes)),
        n_atac_modifier_genes = length(modifier_genes),
        multiome_contributing_genes =
          .rc_ra_collapse(toupper(contribution_genes)),
        n_multiome_contributing_genes = length(contribution_genes),
        has_rna_evidence = length(rna_genes) > 0L,
        has_atac_regulatory_evidence = length(modifier_genes) > 0L,
        has_active_multiome_contribution = length(contribution_genes) > 0L,
        stringsAsFactors = FALSE
      )
    })
    output[[group_index]] <- do.call(rbind, rows)
  }
  answer <- do.call(rbind, output)
  rownames(answer) <- NULL
  answer
}

#' Build reaction formulas, GPR annotations, and evidence provenance
#'
#' Creates one reaction catalog row per GEM reaction and, when Layer 1 is
#' supplied, one evidence-provenance row per condition, cell type, and reaction.
#' The formula is reconstructed directly from the stoichiometric matrix using
#' metabolite names and compartments. Evidence is classified as `RNA+ATAC` only
#' when accessibility-derived regulation changes integrated gene support relative
#' to RNA support in at least one selected unit. Reactions with observed RNA but
#' no active ATAC contribution are labelled `RNA-only`.
#'
#' @param gem A validated RegCompass GEM.
#' @param layer1 Optional Layer 1 result.
#' @param reaction_ids Optional reaction identifiers. The default uses all GEM
#'   reactions.
#' @param condition_col,celltype_col Optional Layer 1 metadata columns.
#' @param evidence_tolerance Numerical tolerance for active evidence.
#' @return A list containing `reactions` and `evidence` data frames.
#' @export
rc_build_reaction_annotations <- function(
    gem, layer1 = NULL, reaction_ids = NULL,
    condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  if (!is.numeric(evidence_tolerance) || length(evidence_tolerance) != 1L ||
      !is.finite(evidence_tolerance) || evidence_tolerance < 0) {
    stop("`evidence_tolerance` must be one non-negative finite number.",
         call. = FALSE)
  }
  catalog <- .rc_ra_reaction_catalog(gem, layer1, reaction_ids)
  evidence <- .rc_ra_group_evidence(
    catalog = catalog,
    layer1 = layer1,
    condition_col = condition_col,
    celltype_col = celltype_col,
    evidence_tolerance = evidence_tolerance
  )
  answer <- list(
    reactions = catalog,
    evidence = evidence,
    params = list(
      condition_col = condition_col,
      celltype_col = celltype_col,
      evidence_tolerance = evidence_tolerance,
      evidence_definition = paste(
        "RNA+ATAC requires gene_support_multiome to differ from",
        "gene_support_rna for at least one GPR gene and unit"
      )
    )
  )
  class(answer) <- c("regcompass_reaction_annotations", "list")
  answer
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

.rc_ra_pairwise_evidence <- function(data, evidence) {
  if (!is.data.frame(data) || !nrow(data) ||
      !is.data.frame(evidence) || !nrow(evidence)) return(data)
  index <- .rc_ra_evidence_index(evidence)
  fields <- c(
    "evidence_class", "rna_supported_genes", "n_rna_supported_genes",
    "atac_modifier_genes", "n_atac_modifier_genes",
    "multiome_contributing_genes", "n_multiome_contributing_genes",
    "has_rna_evidence", "has_atac_regulatory_evidence",
    "has_active_multiome_contribution"
  )
  for (suffix in c("a", "b")) {
    condition <- data[[paste0("condition_", suffix)]]
    key <- paste(data$reaction_id, condition, data$cell_type, sep = "\001")
    idx <- match(key, index)
    for (field in fields) {
      data[[paste0(field, "_", suffix)]] <- evidence[[field]][idx]
    }
  }
  data$evidence_comparison <- paste0(
    data$condition_a, "=", data$evidence_class_a, ";",
    data$condition_b, "=", data$evidence_class_b
  )
  data
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

.rc_ra_omnibus_evidence <- function(data, evidence, conditions) {
  if (!is.data.frame(data) || !nrow(data) ||
      !is.data.frame(evidence) || !nrow(evidence)) return(data)
  index <- .rc_ra_evidence_index(evidence)
  summary <- lapply(seq_len(nrow(data)), function(i) {
    keys <- paste(
      data$reaction_id[[i]], conditions, data$cell_type[[i]], sep = "\001"
    )
    one <- evidence[match(keys, index), , drop = FALSE]
    classes <- as.character(one$evidence_class)
    overall <- if (any(classes == "RNA+ATAC", na.rm = TRUE)) {
      "RNA+ATAC"
    } else if (any(classes == "RNA-only", na.rm = TRUE)) {
      "RNA-only"
    } else if (any(classes == "GPR/no-observed-RNA", na.rm = TRUE)) {
      "GPR/no-observed-RNA"
    } else {
      "structural/no-GPR"
    }
    data.frame(
      evidence_class = overall,
      evidence_by_condition = paste0(
        conditions, "=", ifelse(is.na(classes), "unknown", classes),
        collapse = ";"
      ),
      multiome_supported_conditions = .rc_ra_collapse(
        conditions[classes == "RNA+ATAC"]
      ),
      has_active_multiome_contribution = any(
        one$has_active_multiome_contribution %in% TRUE, na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, summary)
  for (column in colnames(summary)) data[[column]] <- summary[[column]]
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

# Enhanced Stage 6 assembly: calculations remain in the canonical function;
# this layer adds biological reaction metadata and persists the enriched result.
rc_regcompass_step_results <- function(
    grn, metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse")) {
  answer <- .rc_step_results_core(
    grn = grn,
    metacells = metacells,
    meta_modules = meta_modules,
    layer1 = layer1,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    species = species
  )
  params <- metacells$params
  answer <- .rc_ra_attach_to_result(
    result = answer,
    gem = gem,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    answer$reaction_catalog,
    file.path(outdir, "reaction_catalog.tsv.gz")
  )
  .rc_write_tsv_gz(
    answer$reaction_evidence,
    file.path(outdir, "reaction_evidence_by_condition_celltype.tsv.gz")
  )
  saveRDS(answer, file.path(outdir, "regcompass_result.rds"))
  answer
}

# Enhanced statistics output. The tests are delegated unchanged to the canonical
# implementation, then biological reaction annotation is joined by reaction ID.
rc_test_condition_reactions <- function(
    x,
    condition_col = NULL,
    celltype_col = NULL,
    conditions = NULL,
    cell_types = NULL,
    comparisons = NULL,
    reaction_ids = NULL,
    target_directions = NULL,
    medium_scenarios = NULL,
    min_units = 5L,
    include_omnibus = TRUE,
    p_adjust_method = "BH",
    p_adjust_scope = c(
      "celltype_contrast_medium", "celltype_contrast", "celltype", "global"
    ),
    wilcox_correct = FALSE,
    eps = 1e-8,
    vmax_tolerance = 1e-6,
    include_scores = FALSE,
    outdir = NULL) {
  source <- x
  answer <- .rc_condition_statistics_core(
    x = x,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    cell_types = cell_types,
    comparisons = comparisons,
    reaction_ids = reaction_ids,
    target_directions = target_directions,
    medium_scenarios = medium_scenarios,
    min_units = min_units,
    include_omnibus = include_omnibus,
    p_adjust_method = p_adjust_method,
    p_adjust_scope = p_adjust_scope,
    wilcox_correct = wilcox_correct,
    eps = eps,
    vmax_tolerance = vmax_tolerance,
    include_scores = include_scores,
    outdir = outdir
  )
  annotation <- .rc_ra_annotation_from_object(source)
  answer <- .rc_ra_enrich_statistics(answer, annotation)
  if (!is.null(outdir) && !is.null(annotation)) {
    .rc_write_tsv_gz(
      answer$pairwise,
      file.path(outdir, "condition_reaction_pairwise.tsv.gz")
    )
    if (nrow(answer$omnibus)) {
      .rc_write_tsv_gz(
        answer$omnibus,
        file.path(outdir, "condition_reaction_omnibus.tsv.gz")
      )
    }
    .rc_write_tsv_gz(
      answer$reaction_catalog,
      file.path(outdir, "condition_reaction_catalog.tsv.gz")
    )
    .rc_write_tsv_gz(
      answer$reaction_evidence,
      file.path(outdir, "condition_reaction_evidence.tsv.gz")
    )
    saveRDS(answer, file.path(outdir, "condition_reaction_statistics.rds"))
  }
  answer
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

# Add biological labels to the existing single-reaction plot without changing
# its statistical or geometric implementation.
rc_plot_condition_reaction <- function(
    x,
    reaction_id,
    cell_type,
    target_direction = NULL,
    medium_scenario = NULL,
    condition_col = NULL,
    celltype_col = NULL,
    conditions = NULL,
    comparisons = NULL,
    min_units = 5L,
    p_adjust_method = "BH",
    p_adjust_scope = c(
      "celltype_contrast_medium", "celltype_contrast", "celltype", "global"
    ),
    annotation_p = c("p_adj", "p_value"),
    significance_threshold = 0.05,
    show_nonsignificant = FALSE,
    show_omnibus = TRUE,
    point_size = 1.8,
    point_alpha = 0.75,
    jitter_width = 0.12,
    box_width = 0.55,
    bracket_step = 0.12,
    title = NULL,
    y_label = "Reaction support score") {
  plot <- .rc_condition_plot_core(
    x = x,
    reaction_id = reaction_id,
    cell_type = cell_type,
    target_direction = target_direction,
    medium_scenario = medium_scenario,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    comparisons = comparisons,
    min_units = min_units,
    p_adjust_method = p_adjust_method,
    p_adjust_scope = p_adjust_scope,
    annotation_p = annotation_p,
    significance_threshold = significance_threshold,
    show_nonsignificant = show_nonsignificant,
    show_omnibus = show_omnibus,
    point_size = point_size,
    point_alpha = point_alpha,
    jitter_width = jitter_width,
    box_width = box_width,
    bracket_step = bracket_step,
    title = title,
    y_label = y_label
  )
  statistics <- attr(plot, "condition_statistics")
  annotation_row <- data.frame()
  evidence_text <- NULL
  if (is.list(statistics) && is.data.frame(statistics$pairwise)) {
    annotation_row <- statistics$pairwise[
      statistics$pairwise$reaction_id == reaction_id &
        statistics$pairwise$cell_type == cell_type &
        (is.null(target_direction) |
           statistics$pairwise$target_direction == target_direction) &
        (is.null(medium_scenario) |
           statistics$pairwise$medium_scenario == medium_scenario),
      , drop = FALSE
    ]
    if (nrow(annotation_row)) {
      annotation_row <- annotation_row[1L, , drop = FALSE]
      evidence_text <- annotation_row$evidence_comparison[[1L]] %||% NULL
      if (is.null(title) && .rc_ra_nonempty(annotation_row$reaction_name[[1L]])) {
        plot <- plot + ggplot2::labs(
          title = paste0(
            annotation_row$reaction_name[[1L]], " (", reaction_id, ", ",
            annotation_row$target_direction[[1L]], ") in ", cell_type
          )
        )
      }
      plot <- plot + ggplot2::labs(
        caption = .rc_ra_plot_caption(annotation_row, evidence_text)
      ) + ggplot2::theme(
        plot.caption = ggplot2::element_text(hjust = 0, size = 8)
      )
    }
  }
  attr(plot, "reaction_annotation") <- annotation_row
  plot
}

.rc_ra_plot_one_precomputed <- function(
    x, statistics, row_id, cell_type, condition_col, celltype_col,
    conditions, annotation_p, significance_threshold,
    show_nonsignificant, show_omnibus, point_size, point_alpha,
    jitter_width, box_width, bracket_step) {
  microcompass <- .rc_condition_stats_microcompass(x)
  meta <- microcompass$unit_meta
  condition_col <- .rc_condition_stats_column(
    meta, condition_col,
    c("condition", "dataset", "Group", "group", "treatment"),
    "condition_col"
  )
  celltype_col <- .rc_condition_stats_column(
    meta, celltype_col,
    c("cell_type", "celltype", "epithelial_or_stem", "CellType"),
    "celltype_col"
  )
  unit_ids <- .rc_condition_stats_unit_ids(meta)
  meta$.unit_id <- unit_ids
  score <- statistics$score[row_id, , drop = TRUE]
  plot_meta <- meta[match(names(score), meta$.unit_id), , drop = FALSE]
  plot_data <- data.frame(
    unit_id = names(score),
    condition = factor(as.character(plot_meta[[condition_col]]), levels = conditions),
    cell_type = as.character(plot_meta[[celltype_col]]),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  plot_data <- plot_data[
    plot_data$cell_type == cell_type & is.finite(plot_data$score), , drop = FALSE
  ]
  pairwise <- statistics$pairwise[
    statistics$pairwise$row_id == row_id &
      statistics$pairwise$cell_type == cell_type, , drop = FALSE
  ]
  score_min <- min(plot_data$score)
  score_max <- max(plot_data$score)
  annotation_data <- .rc_plot_condition_annotations(
    pairwise = pairwise,
    condition_levels = conditions,
    p_column = annotation_p,
    significance_threshold = significance_threshold,
    show_nonsignificant = show_nonsignificant,
    score_min = score_min,
    score_max = score_max,
    bracket_step = bracket_step
  )
  target <- pairwise[1L, , drop = FALSE]
  omnibus_subtitle <- NULL
  if (isTRUE(show_omnibus) && nrow(statistics$omnibus)) {
    omnibus <- statistics$omnibus[
      statistics$omnibus$row_id == row_id &
        statistics$omnibus$cell_type == cell_type, , drop = FALSE
    ]
    if (nrow(omnibus) == 1L && is.finite(omnibus$p_adj[[1L]])) {
      omnibus_subtitle <- paste0(
        "Kruskal-Wallis ", statistics$params$p_adjust_method,
        "-adjusted P = ",
        format.pval(omnibus$p_adj[[1L]], digits = 3, eps = 1e-4)
      )
    }
  }
  plot <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = condition, y = score, fill = condition)
  ) +
    ggplot2::geom_boxplot(
      width = box_width, outlier.shape = NA, alpha = 0.65, linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = condition), width = jitter_width, height = 0,
      size = point_size, alpha = point_alpha, show.legend = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        target$reaction_name[[1L]], " (", target$reaction_id[[1L]], ", ",
        target$target_direction[[1L]], ") in ", cell_type
      ),
      subtitle = omnibus_subtitle,
      caption = .rc_ra_plot_caption(
        target, target$evidence_comparison[[1L]] %||% NULL
      ),
      x = NULL, y = "Reaction support score", fill = "Condition"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(8, 18, 8, 8),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  if (nrow(annotation_data)) {
    for (i in seq_len(nrow(annotation_data))) {
      one <- annotation_data[i, , drop = FALSE]
      plot <- plot +
        ggplot2::annotate(
          "segment", x = one$xmin, xend = one$xmax,
          y = one$y, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment", x = one$xmin, xend = one$xmin,
          y = one$y - one$tip, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "segment", x = one$xmax, xend = one$xmax,
          y = one$y - one$tip, yend = one$y, linewidth = 0.45
        ) +
        ggplot2::annotate(
          "text", x = (one$xmin + one$xmax) / 2,
          y = one$text_y, label = one$label, vjust = 0, size = 4
        )
    }
    score_range <- score_max - score_min
    if (!is.finite(score_range) || score_range <= 0) {
      score_range <- max(1, abs(score_min), abs(score_max))
    }
    plot <- plot + ggplot2::coord_cartesian(
      ylim = c(
        score_min - score_range * 0.06,
        max(annotation_data$text_y) + score_range * 0.08
      ),
      clip = "off"
    )
  }
  attr(plot, "plot_data") <- plot_data
  attr(plot, "annotation_data") <- annotation_data
  attr(plot, "reaction_annotation") <- target
  plot
}

#' Plot significant condition responses for reactions selected by genes
#'
#' Runs the condition statistics once, selects scored reactions containing the
#' requested metabolic genes, ranks significant reaction-direction targets, and
#' returns a named collection of annotated boxplots.
#'
#' @param x Annotated RegCompass result.
#' @param genes Metabolic gene symbols used to select GPR reactions.
#' @param cell_type One cell type.
#' @param condition_col,celltype_col Metadata columns.
#' @param conditions Ordered conditions.
#' @param comparisons Optional condition pairs.
#' @param target_directions Optional scored directions.
#' @param medium_scenario Optional medium.
#' @param evidence_class Optional group evidence classes used for gene-reaction
#'   selection.
#' @param p_adj_max Maximum adjusted pairwise P value.
#' @param min_abs_rank_biserial Minimum absolute rank-biserial effect.
#' @param max_reactions Maximum number of reaction-direction plots.
#' @param annotation_p Significance label column.
#' @param outdir Optional directory for PDF plots and selection tables.
#' @return A `regcompass_gene_reaction_plots` list.
#' @export
rc_plot_condition_gene_reactions <- function(
    x, genes, cell_type,
    condition_col = NULL, celltype_col = NULL,
    conditions = NULL, comparisons = NULL,
    target_directions = NULL, medium_scenario = NULL,
    evidence_class = NULL,
    p_adj_max = 0.05, min_abs_rank_biserial = 0.30,
    max_reactions = 12L,
    annotation_p = c("p_adj", "p_value"),
    significance_threshold = 0.05,
    show_nonsignificant = FALSE, show_omnibus = TRUE,
    point_size = 1.8, point_alpha = 0.75,
    jitter_width = 0.12, box_width = 0.55, bracket_step = 0.12,
    outdir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required.", call. = FALSE)
  }
  annotation_p <- match.arg(annotation_p)
  if (!is.numeric(p_adj_max) || length(p_adj_max) != 1L ||
      !is.finite(p_adj_max) || p_adj_max <= 0 || p_adj_max > 1 ||
      !is.numeric(min_abs_rank_biserial) ||
      length(min_abs_rank_biserial) != 1L ||
      !is.finite(min_abs_rank_biserial) || min_abs_rank_biserial < 0 ||
      min_abs_rank_biserial > 1 ||
      !is.numeric(max_reactions) || length(max_reactions) != 1L ||
      !is.finite(max_reactions) || max_reactions < 1) {
    stop("Gene-reaction plot selection thresholds are invalid.", call. = FALSE)
  }
  max_reactions <- as.integer(max_reactions)
  selection <- rc_select_gene_reactions(
    x = x,
    genes = genes,
    cell_types = cell_type,
    evidence_class = evidence_class
  )
  if (!length(selection$reaction_ids)) {
    stop("No scored reactions matched the requested genes and evidence filters.",
         call. = FALSE)
  }
  statistics <- rc_test_condition_reactions(
    x = x,
    condition_col = condition_col,
    celltype_col = celltype_col,
    conditions = conditions,
    cell_types = cell_type,
    comparisons = comparisons,
    min_units = 5L,
    include_omnibus = TRUE,
    p_adjust_method = "BH",
    p_adjust_scope = "celltype_contrast_medium",
    include_scores = TRUE
  )
  pairwise <- statistics$pairwise[
    statistics$pairwise$reaction_id %in% selection$reaction_ids &
      statistics$pairwise$cell_type == cell_type &
      is.finite(statistics$pairwise$p_adj) &
      statistics$pairwise$p_adj <= p_adj_max &
      is.finite(statistics$pairwise$rank_biserial_b_minus_a) &
      abs(statistics$pairwise$rank_biserial_b_minus_a) >=
        min_abs_rank_biserial,
    , drop = FALSE
  ]
  if (!is.null(target_directions)) {
    pairwise <- pairwise[
      pairwise$target_direction %in% as.character(target_directions), , drop = FALSE
    ]
  }
  if (!is.null(medium_scenario)) {
    pairwise <- pairwise[
      pairwise$medium_scenario %in% as.character(medium_scenario), , drop = FALSE
    ]
  }
  if (!nrow(pairwise)) {
    stop("No gene-associated targets passed the significance and effect filters.",
         call. = FALSE)
  }
  split_rows <- split(seq_len(nrow(pairwise)), pairwise$row_id)
  ranked <- do.call(rbind, lapply(split_rows, function(rows) {
    one <- pairwise[rows, , drop = FALSE]
    data.frame(
      row_id = one$row_id[[1L]],
      reaction_id = one$reaction_id[[1L]],
      reaction_name = one$reaction_name[[1L]],
      target_direction = one$target_direction[[1L]],
      medium_scenario = one$medium_scenario[[1L]],
      tested_formula = one$tested_formula[[1L]],
      genes = one$genes[[1L]],
      min_p_adj = min(one$p_adj, na.rm = TRUE),
      max_abs_rank_biserial = max(
        abs(one$rank_biserial_b_minus_a), na.rm = TRUE
      ),
      max_abs_delta_median = max(
        abs(one$delta_median_score_b_minus_a), na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }))
  ranked <- ranked[order(
    ranked$min_p_adj,
    -ranked$max_abs_rank_biserial,
    -ranked$max_abs_delta_median,
    ranked$reaction_id,
    ranked$target_direction
  ), , drop = FALSE]
  ranked <- utils::head(ranked, max_reactions)
  if (is.null(conditions)) conditions <- statistics$params$conditions
  plots <- lapply(ranked$row_id, function(row_id) {
    .rc_ra_plot_one_precomputed(
      x = x,
      statistics = statistics,
      row_id = row_id,
      cell_type = cell_type,
      condition_col = condition_col,
      celltype_col = celltype_col,
      conditions = conditions,
      annotation_p = annotation_p,
      significance_threshold = significance_threshold,
      show_nonsignificant = show_nonsignificant,
      show_omnibus = show_omnibus,
      point_size = point_size,
      point_alpha = point_alpha,
      jitter_width = jitter_width,
      box_width = box_width,
      bracket_step = bracket_step
    )
  })
  names(plots) <- ranked$row_id
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    .rc_write_tsv_gz(
      ranked,
      file.path(outdir, "selected_gene_reaction_targets.tsv.gz")
    )
    .rc_write_tsv_gz(
      pairwise,
      file.path(outdir, "selected_gene_reaction_pairwise.tsv.gz")
    )
    for (i in seq_along(plots)) {
      safe <- gsub("[^A-Za-z0-9_.-]+", "_", ranked$row_id[[i]])
      ggplot2::ggsave(
        filename = file.path(outdir, paste0(safe, ".pdf")),
        plot = plots[[i]], width = 7, height = 5.5
      )
    }
  }
  answer <- list(
    plots = plots,
    selected_targets = ranked,
    pairwise_hits = pairwise,
    statistics = statistics,
    gene_selection = selection
  )
  class(answer) <- c("regcompass_gene_reaction_plots", "list")
  answer
}
