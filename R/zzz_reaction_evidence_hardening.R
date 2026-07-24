# Reaction-level evidence hardening.
# Multiome evidence is assigned only when ATAC integration changes the GPR-
# aggregated reaction capacity, not merely when one participating gene changes.

.rc_ra_reaction_capacity_pair <- function(catalog, layer1) {
  empty <- list(rna = NULL, multiome = NULL, resolution = "unavailable")
  if (!is.list(layer1) || !is.list(layer1$parsed_gpr) ||
      is.null(layer1$gene_support_rna) ||
      is.null(layer1$gene_support_multiome)) return(empty)
  reaction_ids <- intersect(
    as.character(catalog$reaction_id), names(layer1$parsed_gpr)
  )
  if (!length(reaction_ids)) return(empty)
  parsed <- layer1$parsed_gpr[reaction_ids]
  params <- layer1$capacity_params %||% list()
  calculate <- function(gene_support) {
    rc_reaction_capacity(
      gpr_list = parsed,
      gene_score = gene_support,
      promiscuity_mode = params$promiscuity_mode %||% "none",
      tau = params$tau %||% 0.20,
      and_method = params$and_method %||% "boltzmann",
      or_method = params$or_method %||% "sum",
      BPPARAM = FALSE
    )
  }
  rna <- calculate(layer1$gene_support_rna)
  multiome <- layer1$reaction_expression
  valid_multiome <- is.matrix(multiome) &&
    all(reaction_ids %in% rownames(multiome)) &&
    identical(colnames(rna), colnames(multiome))
  if (valid_multiome) {
    multiome <- multiome[reaction_ids, colnames(rna), drop = FALSE]
  } else {
    multiome <- calculate(layer1$gene_support_multiome)
  }
  list(
    rna = rna[reaction_ids, , drop = FALSE],
    multiome = multiome[reaction_ids, , drop = FALSE],
    resolution = "reaction_capacity"
  )
}

.rc_ra_group_evidence <- function(
    catalog, layer1, condition_col = NULL, celltype_col = NULL,
    evidence_tolerance = 1e-8) {
  empty <- data.frame(
    reaction_id = character(), condition = character(), cell_type = character(),
    evidence_class = character(), evidence_resolution = character(),
    n_units = integer(),
    rna_supported_genes = character(), n_rna_supported_genes = integer(),
    atac_modifier_genes = character(), n_atac_modifier_genes = integer(),
    multiome_contributing_genes = character(),
    n_multiome_contributing_genes = integer(),
    has_rna_evidence = logical(),
    has_atac_regulatory_evidence = logical(),
    has_active_multiome_contribution = logical(),
    median_multiome_capacity_shift = numeric(),
    max_abs_multiome_capacity_shift = numeric(),
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
  capacities <- .rc_ra_reaction_capacity_pair(catalog, layer1)
  if (!is.null(capacities$rna)) {
    capacities$rna <- capacities$rna[, colnames(rna), drop = FALSE]
    capacities$multiome <- capacities$multiome[, colnames(rna), drop = FALSE]
  }
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
      all_genes <- gene_lists[[reaction_id]]
      genes <- intersect(all_genes, rownames(rna))
      flag_genes <- function(flag) {
        if (!length(genes)) return(character())
        values <- unname(flag[genes])
        values[is.na(values)] <- FALSE
        genes[values]
      }
      rna_genes <- flag_genes(rna_flag)
      modifier_genes <- flag_genes(modifier_flag)
      contribution_genes <- flag_genes(contribution_flag)
      has_rna_capacity <- FALSE
      has_capacity_shift <- FALSE
      median_shift <- NA_real_
      max_abs_shift <- NA_real_
      evidence_resolution <- capacities$resolution
      if (!is.null(capacities$rna) &&
          reaction_id %in% rownames(capacities$rna)) {
        rna_capacity <- capacities$rna[reaction_id, keep, drop = TRUE]
        multiome_capacity <- capacities$multiome[reaction_id, keep, drop = TRUE]
        has_rna_capacity <- any(
          is.finite(rna_capacity) & rna_capacity > evidence_tolerance
        )
        valid <- is.finite(rna_capacity) & is.finite(multiome_capacity)
        shifts <- multiome_capacity[valid] - rna_capacity[valid]
        if (length(shifts)) {
          has_capacity_shift <- any(abs(shifts) > evidence_tolerance)
          median_shift <- stats::median(shifts)
          max_abs_shift <- max(abs(shifts))
        }
      } else {
        evidence_resolution <- "gene_support_fallback"
        has_rna_capacity <- length(rna_genes) > 0L
        has_capacity_shift <- length(contribution_genes) > 0L
      }
      evidence_class <- if (!length(all_genes)) {
        "structural/no-GPR"
      } else if (!has_rna_capacity) {
        "GPR/no-observed-RNA"
      } else if (has_capacity_shift) {
        "RNA+ATAC"
      } else {
        "RNA-only"
      }
      data.frame(
        reaction_id = reaction_id,
        condition = condition,
        cell_type = cell_type,
        evidence_class = evidence_class,
        evidence_resolution = evidence_resolution,
        n_units = sum(keep),
        rna_supported_genes = .rc_ra_collapse(toupper(rna_genes)),
        n_rna_supported_genes = length(rna_genes),
        atac_modifier_genes = .rc_ra_collapse(toupper(modifier_genes)),
        n_atac_modifier_genes = length(modifier_genes),
        multiome_contributing_genes =
          .rc_ra_collapse(toupper(contribution_genes)),
        n_multiome_contributing_genes = length(contribution_genes),
        has_rna_evidence = has_rna_capacity,
        has_atac_regulatory_evidence = length(modifier_genes) > 0L,
        has_active_multiome_contribution = has_capacity_shift,
        median_multiome_capacity_shift = median_shift,
        max_abs_multiome_capacity_shift = max_abs_shift,
        stringsAsFactors = FALSE
      )
    })
    output[[group_index]] <- do.call(rbind, rows)
  }
  answer <- do.call(rbind, output)
  rownames(answer) <- NULL
  answer
}
