from pathlib import Path
import re

ROOT = Path('.')
R = ROOT / 'R'
TESTS = ROOT / 'tests' / 'testthat'


def function_span(text, name):
    pattern = re.compile(rf'(?m)^{re.escape(name)}\s*<-\s*function\s*\(')
    match = pattern.search(text)
    if match is None:
        raise RuntimeError(f'Function not found: {name}')
    brace = text.find('{', match.end())
    if brace < 0:
        raise RuntimeError(f'Opening brace not found: {name}')
    depth = 0
    quote = None
    escaped = False
    comment = False
    i = brace
    while i < len(text):
        ch = text[i]
        if comment:
            if ch == '\n':
                comment = False
            i += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif ch == '\\':
                escaped = True
            elif ch == quote:
                quote = None
            i += 1
            continue
        if ch in ('"', "'"):
            quote = ch
        elif ch == '#':
            comment = True
        elif ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = i + 1
                while end < len(text) and text[end] in ' \t':
                    end += 1
                while end < len(text) and text[end] == '\n':
                    end += 1
                return match.start(), end
        i += 1
    raise RuntimeError(f'Function did not close: {name}')


def replace_function(path, name, replacement):
    text = path.read_text()
    start, end = function_span(text, name)
    path.write_text(text[:start] + replacement.rstrip() + '\n\n' + text[end:])


def replace_once(path, old, new):
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f'{path}: expected one match, found {count}: {old[:100]}')
    path.write_text(text.replace(old, new, 1))


layer2_direction = r'''.rc_layer2_direction_contract <- function(x) {
  tab <- as.data.frame(x$lp_diagnostics)
  required <- c(
    "reaction_id", "target_direction", "medium_scenario", "row_id"
  )
  if (identical(as.character(x$model_mode), "meta_module_gem")) {
    required <- c("cell_type", required)
  }
  if (!all(required %in% colnames(tab))) {
    stop("Layer 2 direction diagnostics are incomplete.", call. = FALSE)
  }
  tab <- unique(tab[, required, drop = FALSE])
  order_cols <- intersect(
    c("cell_type", "medium_scenario", "reaction_id",
      "target_direction", "row_id"),
    colnames(tab)
  )
  tab[do.call(order, tab[order_cols]), , drop = FALSE]
}
'''
replace_function(R / 'step_layer2.R', '.rc_layer2_direction_contract', layer2_direction)

comparison = r'''.rc_layer2_comparison_table <- function(
    layer2, layer1, condition_col, celltype_col) {
  penalty <- layer2$penalty
  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  unit_meta <- layer2$unit_meta
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else {
    as.character(unit_meta$pool_id)
  }
  unit_meta <- unit_meta[match(colnames(penalty), unit_id), , drop = FALSE]
  unit_celltype <- as.character(unit_meta[[celltype_col]])
  grid <- expand.grid(
    row_index = seq_len(nrow(penalty)),
    unit_index = seq_len(ncol(penalty)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scoped_celltype <- row_meta$cell_type[grid$row_index]
  matching_scope <- is.na(scoped_celltype) |
    scoped_celltype == unit_celltype[grid$unit_index]
  grid <- grid[matching_scope, , drop = FALSE]
  if (!nrow(grid)) {
    stop("No Layer 2 rows match their cell-type units.", call. = FALSE)
  }
  reaction <- row_meta$reaction_id[grid$row_index]
  unit <- colnames(penalty)[grid$unit_index]
  index_matrix <- function(x) {
    as.numeric(x[cbind(grid$row_index, grid$unit_index)])
  }
  reaction_unit_matrix <- function(x) {
    if (!is.numeric(x) || is.null(dim(x)) ||
        !all(unique(reaction) %in% rownames(x)) ||
        !identical(colnames(x), colnames(penalty))) {
      stop("Layer 1 reaction diagnostics are not aligned.", call. = FALSE)
    }
    as.numeric(x[cbind(match(reaction, rownames(x)), grid$unit_index)])
  }
  omega <- layer2$params$omega
  vmax <- index_matrix(layer2$vmax)
  primary <- index_matrix(layer2$penalty_condition_full_oof)
  normalized <- primary / (omega * vmax)
  normalized[!is.finite(primary) | !is.finite(vmax) | vmax <= 0] <- NA_real_
  data.frame(
    row_id = rownames(penalty)[grid$row_index],
    reaction_id = reaction,
    direction = row_meta$target_direction[grid$row_index],
    medium = row_meta$medium_scenario[grid$row_index],
    cell_type = unit_celltype[grid$unit_index],
    condition = as.character(
      unit_meta[[condition_col]][grid$unit_index]
    ),
    metacell_id = unit,
    penalty_condition_full_oof = primary,
    penalty_common_oof = index_matrix(layer2$penalty_common_oof),
    penalty_condition_unique_increment = index_matrix(
      layer2$penalty_condition_unique_increment
    ),
    penalty_rna_only = index_matrix(layer2$penalty_rna_only),
    penalty_per_target_flux = normalized,
    vmax = vmax,
    condition_full_oof_available = is.finite(primary),
    condition_full_support_fraction = reaction_unit_matrix(
      layer1$reaction_condition_full_support_fraction
    ),
    common_support_fraction = reaction_unit_matrix(
      layer1$reaction_common_support_fraction
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    comparability_class =
      "same_celltype_conditions_on_one_celltype_medium_union_gem",
    stringsAsFactors = FALSE
  )
}
'''
replace_function(R / 'step_layer2.R', '.rc_layer2_comparison_table', comparison)

condition_comparison = r'''.rc_condition_penalty_comparison <- function(
    microcompass, condition_col = "condition", celltype_col = "cell_type",
    eps = 1e-8, vmax_tolerance = 1e-6) {
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(vmax_tolerance) || length(vmax_tolerance) != 1L ||
      !is.finite(vmax_tolerance) || vmax_tolerance < 0) {
    stop("`vmax_tolerance` must be one finite non-negative number.",
         call. = FALSE)
  }
  penalty <- as.matrix(microcompass$penalty)
  vmax <- as.matrix(microcompass$vmax)
  meta <- microcompass$unit_meta
  required <- c(condition_col, celltype_col)
  valid_matrix <- function(x) {
    is.numeric(x) && !is.null(rownames(x)) && !is.null(colnames(x)) &&
      !anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x))
  }
  if (!valid_matrix(penalty) || !valid_matrix(vmax) ||
      !setequal(rownames(penalty), rownames(vmax)) ||
      !setequal(colnames(penalty), colnames(vmax))) {
    stop("microCOMPASS penalty and vmax matrices must align.", call. = FALSE)
  }
  vmax <- vmax[rownames(penalty), colnames(penalty), drop = FALSE]
  if (!is.data.frame(meta) || !all(required %in% colnames(meta))) {
    stop("microCOMPASS unit metadata lack condition/cell-type columns.",
         call. = FALSE)
  }
  unit_id <- if ("unit_id" %in% colnames(meta)) {
    as.character(meta$unit_id)
  } else if ("pool_id" %in% colnames(meta)) {
    as.character(meta$pool_id)
  } else {
    stop("microCOMPASS unit metadata lack unit_id/pool_id.", call. = FALSE)
  }
  if (anyNA(unit_id) || any(!nzchar(trimws(unit_id))) || anyDuplicated(unit_id) ||
      !setequal(colnames(penalty), unit_id)) {
    stop("microCOMPASS units and metadata are invalid or different.",
         call. = FALSE)
  }
  meta$unit_id <- unit_id
  meta <- meta[match(colnames(penalty), meta$unit_id), , drop = FALSE]
  condition_value <- trimws(as.character(meta[[condition_col]]))
  celltype_value <- trimws(as.character(meta[[celltype_col]]))
  if (anyNA(condition_value) || any(!nzchar(condition_value)) ||
      anyNA(celltype_value) || any(!nzchar(celltype_value))) {
    stop("microCOMPASS condition/cell-type metadata are incomplete.",
         call. = FALSE)
  }

  omega <- microcompass$params$omega %||% 0.95
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("microCOMPASS `omega` must be in (0, 1].", call. = FALSE)
  }
  vmax_invariant <- vapply(seq_len(nrow(vmax)), function(i) {
    values <- vmax[i, is.finite(vmax[i, ]), drop = TRUE]
    if (length(values) <= 1L) return(TRUE)
    diff(range(values)) <= vmax_tolerance *
      max(1, abs(stats::median(values)))
  }, logical(1))
  if (any(!vmax_invariant)) {
    stop(
      "Target vmax differs among units assigned to the same structural row: ",
      paste(utils::head(rownames(vmax)[!vmax_invariant], 10L), collapse = ", "),
      call. = FALSE
    )
  }
  required_target_flux <- omega * vmax
  normalized <- matrix(
    NA_real_, nrow(penalty), ncol(penalty), dimnames = dimnames(penalty)
  )
  valid <- is.finite(penalty) & is.finite(required_target_flux) &
    required_target_flux > 0
  normalized[valid] <- penalty[valid] / required_target_flux[valid]

  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  row_meta$row_id <- rownames(penalty)
  strata <- unique(meta[, c(condition_col, celltype_col), drop = FALSE])
  summary_rows <- lapply(seq_len(nrow(strata)), function(i) {
    condition <- as.character(strata[[condition_col]][[i]])
    cell_type <- as.character(strata[[celltype_col]][[i]])
    column_keep <- condition_value == condition & celltype_value == cell_type
    row_keep <- is.na(row_meta$cell_type) | row_meta$cell_type == cell_type
    if (!any(row_keep)) return(NULL)
    row_index <- which(row_keep)
    median_matrix <- function(x) {
      value <- matrixStats::rowMedians(
        x[row_index, column_keep, drop = FALSE], na.rm = TRUE
      )
      value[is.nan(value)] <- NA_real_
      value
    }
    median_penalty <- median_matrix(penalty)
    median_vmax <- median_matrix(vmax)
    median_required_flux <- median_matrix(required_target_flux)
    median_normalized <- median_matrix(normalized)
    data.frame(
      row_id = row_meta$row_id[row_index],
      reaction_id = row_meta$reaction_id[row_index],
      target_direction = row_meta$target_direction[row_index],
      medium_scenario = row_meta$medium_scenario[row_index],
      condition = condition,
      cell_type = cell_type,
      median_penalty = median_penalty,
      median_vmax = median_vmax,
      median_required_target_flux = median_required_flux,
      median_penalty_per_target_flux = median_normalized,
      support_score = -log(median_normalized + eps),
      priority_rank = NA_integer_,
      ranking_metric = "minimum_penalty_per_required_target_flux",
      ranking_scope = "condition_x_celltype_x_medium",
      n_metacells = sum(column_keep),
      stringsAsFactors = FALSE
    )
  })
  summary_rows <- Filter(Negate(is.null), summary_rows)
  if (!length(summary_rows)) {
    stop("No condition rows match their cell-type structural models.",
         call. = FALSE)
  }
  ranking <- do.call(rbind, summary_rows)
  rank_group <- interaction(
    ranking$cell_type, ranking$condition, ranking$medium_scenario,
    drop = TRUE, lex.order = TRUE
  )
  for (rows in split(seq_len(nrow(ranking)), rank_group)) {
    ranking$priority_rank[rows] <- as.integer(rank(
      ranking$median_penalty_per_target_flux[rows],
      ties.method = "min", na.last = "keep"
    ))
  }
  ranking <- ranking[order(
    ranking$cell_type, ranking$condition, ranking$medium_scenario,
    ranking$priority_rank, ranking$reaction_id,
    ranking$target_direction, na.last = TRUE
  ), , drop = FALSE]
  rownames(ranking) <- NULL

  contrast_rows <- list()
  index <- 0L
  for (cell_type in unique(ranking$cell_type)) {
    one <- ranking[ranking$cell_type == cell_type, , drop = FALSE]
    conditions <- unique(as.character(one$condition))
    if (length(conditions) < 2L) next
    for (pair in utils::combn(conditions, 2L, simplify = FALSE)) {
      a <- one[one$condition == pair[[1L]], , drop = FALSE]
      b <- one[one$condition == pair[[2L]], , drop = FALSE]
      b <- b[match(a$row_id, b$row_id), , drop = FALSE]
      if (anyNA(b$row_id)) {
        stop("Conditions within a cell type contain different reaction rows.",
             call. = FALSE)
      }
      index <- index + 1L
      contrast_rows[[index]] <- data.frame(
        row_id = a$row_id,
        reaction_id = a$reaction_id,
        target_direction = a$target_direction,
        medium_scenario = a$medium_scenario,
        cell_type = cell_type,
        condition_a = pair[[1L]],
        condition_b = pair[[2L]],
        median_penalty_a = a$median_penalty,
        median_penalty_b = b$median_penalty,
        median_penalty_per_target_flux_a =
          a$median_penalty_per_target_flux,
        median_penalty_per_target_flux_b =
          b$median_penalty_per_target_flux,
        priority_rank_a = a$priority_rank,
        priority_rank_b = b$priority_rank,
        delta_support_b_minus_a = b$support_score - a$support_score,
        higher_supported_condition = ifelse(
          b$support_score > a$support_score, pair[[2L]],
          ifelse(a$support_score > b$support_score, pair[[1L]], "tie")
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  contrast <- if (length(contrast_rows)) do.call(rbind, contrast_rows) else
    data.frame()
  if (nrow(contrast)) rownames(contrast) <- NULL
  list(
    summary = ranking,
    ranking = ranking,
    contrast = contrast,
    analysis_mode = if (length(unique(ranking$condition)) == 1L) {
      "single_condition_reaction_ranking"
    } else {
      "multi_condition_reaction_ranking_and_pairwise_comparison"
    },
    ranking_formula = "penalty / (omega * vmax)",
    ranking_scope = "condition_x_celltype_x_medium",
    structural_scope = "cell_type_x_medium",
    comparison_workflow = paste(
      "conditions are compared only within one cell type on that cell type's",
      "medium-specific union GEM"
    )
  )
}
'''
replace_function(R / 'condition_layer1.R', '.rc_condition_penalty_comparison', condition_comparison)

# Condition statistics must only test rows belonging to the selected cell type.
stats_path = R / 'condition_statistics.R'
replace_once(
    stats_path,
    '  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))\n  row_meta$row_id <- rownames(penalty)\n  row_keep <- rep(TRUE, nrow(row_meta))\n',
    '  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))\n'
    '  row_meta$row_id <- rownames(penalty)\n'
    '  scoped_types <- unique(row_meta$cell_type[!is.na(row_meta$cell_type)])\n'
    '  unknown_scoped <- setdiff(scoped_types, available_cell_types)\n'
    '  if (length(unknown_scoped)) {\n'
    '    stop("Reaction rows contain unknown cell types: ",\n'
    '         paste(unknown_scoped, collapse = ", "), call. = FALSE)\n'
    '  }\n'
    '  row_keep <- is.na(row_meta$cell_type) |\n'
    '    row_meta$cell_type %in% cell_types\n'
)
replace_once(
    stats_path,
    '      for (i in seq_len(nrow(score))) {\n',
    '      row_indices <- which(is.na(row_meta$cell_type) |\n'
    '        row_meta$cell_type == cell_type)\n'
    '      for (i in row_indices) {\n'
)
# Replace the second loop occurrence for omnibus.
text = stats_path.read_text()
old = '      for (i in seq_len(nrow(score))) {\n'
pos = text.find(old)
if pos < 0:
    raise RuntimeError('condition statistics omnibus loop not found')
text = text[:pos] + (
    '      row_indices <- which(is.na(row_meta$cell_type) |\n'
    '        row_meta$cell_type == cell_type)\n'
    '      for (i in row_indices) {\n'
) + text[pos + len(old):]
text = text.replace(
    'under one shared structural GEM.',
    'under one structural GEM specific to that cell type and medium.'
)
stats_path.write_text(text)

# A reaction-ID anchor may be absent from one cell type; do not let that prevent
# valid anchors from being evaluated in other cell types.
target_path = R / 'target_union.R'
replace_once(
    target_path,
    '    if (!length(gene_reactions)) {\n      stop(\n        "The selected genes do not resolve to original Layer 2 core targets.",\n        call. = FALSE\n      )\n    }\n',
    '    if (!length(gene_reactions) && !length(requested_reactions)) {\n'
    '      stop(\n'
    '        "The selected genes do not resolve to original Layer 2 core targets.",\n'
    '        call. = FALSE\n'
    '      )\n'
    '    }\n'
)
replace_once(
    target_path,
    '  reactions <- union(requested_reactions, gene_reactions)\n',
    '  reactions <- union(requested_reactions, gene_reactions)\n'
    '  if (!length(reactions)) {\n'
    '    stop("No selected anchor is available in this cell type.",\n'
    '         call. = FALSE)\n'
    '  }\n'
)

target_definition = r'''.rc_build_target_union_definition <- function(
    gem, merged_core_reactions, merged_reaction_membership,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    cached_reaction_ids = NULL,
    celltype_col = "cell_type") {
  gene_match <- match.arg(gene_match)
  required <- c(celltype_col, "reaction_id")
  if (!is.data.frame(merged_core_reactions) ||
      !all(required %in% colnames(merged_core_reactions)) ||
      !is.data.frame(merged_reaction_membership) ||
      !all(required %in% colnames(merged_reaction_membership))) {
    stop("Merged meta-module tables must be scoped by cell type.",
         call. = FALSE)
  }
  if (!is.list(cached_reaction_ids) || is.null(names(cached_reaction_ids))) {
    stop("Cached reaction IDs must be a named list by cell type.",
         call. = FALSE)
  }
  core_types <- sort(unique(trimws(as.character(
    merged_core_reactions[[celltype_col]]
  ))))
  membership_types <- sort(unique(trimws(as.character(
    merged_reaction_membership[[celltype_col]]
  ))))
  if (!setequal(core_types, membership_types) ||
      !setequal(core_types, names(cached_reaction_ids))) {
    stop("Core, membership and cached GEMs cover different cell types.",
         call. = FALSE)
  }

  build_one <- function(cell_type) {
    cores <- merged_core_reactions[
      trimws(as.character(merged_core_reactions[[celltype_col]])) == cell_type,
      , drop = FALSE
    ]
    membership <- merged_reaction_membership[
      trimws(as.character(merged_reaction_membership[[celltype_col]])) ==
        cell_type, , drop = FALSE
    ]
    available <- .rc_target_union_normalize_ids(
      cached_reaction_ids[[cell_type]]
    )
    requested <- .rc_target_union_normalize_ids(core_reaction_ids)
    requested <- intersect(requested, union(cores$reaction_id, available))
    selected <- tryCatch(
      .rc_target_union_core_rows(
        gem = gem,
        available_core_reactions = cores$reaction_id,
        core_reaction_ids = requested,
        core_genes = core_genes,
        gene_match = gene_match
      ),
      error = function(e) NULL
    )
    if (is.null(selected) || !nrow(selected)) return(NULL)
    selected <- selected[selected$reaction_id %in% available, , drop = FALSE]
    if (!nrow(selected)) return(NULL)
    selected[[celltype_col]] <- cell_type
    selected$catalogue_id <- paste0(
      "CELLTYPE_META_MODULES::",
      utils::URLencode(cell_type, reserved = TRUE)
    )
    selected_core <- selected[selected$is_core %in% TRUE, , drop = FALSE]
    selected_noncore <- selected[!selected$is_core %in% TRUE, , drop = FALSE]

    catalogue <- .rc_target_union_direct_crossref_relations(gem, selected)
    if (!nrow(catalogue)) return(NULL)
    catalogue[[celltype_col]] <- cell_type
    catalogue$catalogue_id <- selected$catalogue_id[[1L]]
    catalogue$anchor_reaction_id <- catalogue$anchor_core_reaction_id
    anchor_core <- stats::setNames(
      as.logical(selected$is_core), as.character(selected$reaction_id)
    )
    catalogue$anchor_is_original_core <- unname(
      anchor_core[as.character(catalogue$anchor_reaction_id)]
    )
    catalogue$anchor_is_original_core[
      is.na(catalogue$anchor_is_original_core)
    ] <- FALSE

    membership_match <- match(
      catalogue$reaction_id, membership$reaction_id
    )
    catalogue$available_in_all_cached_union_gems <-
      catalogue$reaction_id %in% available
    catalogue$present_in_merged_catalogue <- !is.na(membership_match)
    catalogue$merged_catalogue_is_core <-
      catalogue$reaction_id %in% cores$reaction_id
    catalogue$merged_catalogue_inclusion_stage <- if (
      "inclusion_stage" %in% colnames(membership)
    ) {
      as.character(membership$inclusion_stage[membership_match])
    } else {
      rep(NA_character_, nrow(catalogue))
    }
    support_only <- catalogue$available_in_all_cached_union_gems &
      !catalogue$present_in_merged_catalogue
    catalogue$merged_catalogue_inclusion_stage[support_only] <-
      "celltype_fastcore_support_not_in_biological_catalogue"
    catalogue$score_target <- !catalogue$merged_catalogue_is_core &
      catalogue$available_in_all_cached_union_gems
    catalogue$target_role <- ifelse(
      catalogue$merged_catalogue_is_core,
      "celltype_core_not_rescored",
      ifelse(
        catalogue$available_in_all_cached_union_gems,
        "direct_database_crossref_noncore",
        "absent_from_celltype_cached_union_gem"
      )
    )
    catalogue$lp_exclusion_reason <- ifelse(
      catalogue$merged_catalogue_is_core,
      "already_scored_in_original_celltype_layer2",
      ifelse(
        catalogue$available_in_all_cached_union_gems,
        NA_character_,
        "absent_from_one_or_more_celltype_medium_union_gems"
      )
    )
    target_relations <- catalogue[catalogue$score_target, , drop = FALSE]
    targets <- if (nrow(target_relations)) {
      .rc_target_union_aggregate_targets(
        target_relations, celltype_col = celltype_col
      )
    } else {
      data.frame()
    }
    if (nrow(targets)) {
      targets$anchor_reaction_ids <- targets$anchor_core_reaction_ids
    }
    expansion_policy <- if (nrow(selected_noncore)) {
      "direct_from_celltype_anchors_via_kegg_reactome_master_rhea_only"
    } else {
      "direct_from_celltype_core_via_kegg_reactome_master_rhea_only"
    }
    summary <- data.frame(
      cell_type = cell_type,
      n_selected_anchors = nrow(selected),
      n_selected_core = nrow(selected_core),
      n_selected_noncore_anchors = nrow(selected_noncore),
      n_direct_crossref_relations = nrow(catalogue),
      n_direct_crossref_reactions = length(unique(catalogue$reaction_id)),
      n_celltype_core_reactions_not_rescored = length(unique(
        catalogue$reaction_id[catalogue$merged_catalogue_is_core]
      )),
      n_cached_union_unavailable_reactions = length(unique(
        catalogue$reaction_id[!catalogue$available_in_all_cached_union_gems]
      )),
      n_expanded_score_targets = nrow(targets),
      n_celltype_catalogue_reactions =
        length(unique(membership$reaction_id)),
      n_reactions_shared_by_celltype_medium_union_gems = length(available),
      gene_match = gene_match,
      expansion_policy = expansion_policy,
      scoring_policy =
        "celltype_noncore_reactions_present_in_all_media_for_that_celltype",
      model_policy =
        "reuse_exact_celltype_medium_union_gems_without_fastcore_rerun",
      stringsAsFactors = FALSE
    )
    list(
      selected_anchor_reactions = selected,
      selected_core_reactions = selected_core,
      selected_noncore_reactions = selected_noncore,
      expanded_reaction_catalog = catalogue,
      expanded_scoring_targets = targets,
      merged_catalogue_membership = membership,
      summary = summary,
      params = list(
        cell_type = cell_type,
        score_targets = if (nrow(targets)) {
          targets[, c(celltype_col, "reaction_id"), drop = FALSE]
        } else {
          data.frame()
        },
        expansion_policy = expansion_policy
      )
    )
  }

  pieces <- Filter(Negate(is.null), lapply(core_types, build_one))
  if (!length(pieces)) {
    stop("No targeted reaction can be mapped within any cell type.",
         call. = FALSE)
  }
  bind <- function(field) .rc_bind_frames_fill(lapply(pieces, `[[`, field))
  targets <- bind("expanded_scoring_targets")
  if (!nrow(targets)) {
    stop(
      "No directly linked non-core target is present in its cell-type union GEMs.",
      call. = FALSE
    )
  }
  selected <- bind("selected_anchor_reactions")
  catalogue <- bind("expanded_reaction_catalog")
  list(
    selected_anchor_reactions = selected,
    selected_core_reactions = bind("selected_core_reactions"),
    selected_noncore_reactions = bind("selected_noncore_reactions"),
    expanded_reaction_catalog = catalogue,
    expanded_scoring_targets = targets,
    merged_catalogue_membership = bind("merged_catalogue_membership"),
    summary = bind("summary"),
    params = list(
      gene_match = gene_match,
      celltype_col = celltype_col,
      cell_types = sort(unique(as.character(selected[[celltype_col]]))),
      score_targets = unique(
        targets[, c(celltype_col, "reaction_id"), drop = FALSE]
      ),
      expansion_policy =
        "direct_database_crossrefs_evaluated_within_cell_type_only"
    )
  )
}
'''
replace_function(target_path, '.rc_build_target_union_definition', target_definition)

aggregate_targets = r'''.rc_target_union_aggregate_targets <- function(
    catalogue, celltype_col = "cell_type") {
  if (!all(c(celltype_col, "reaction_id") %in% colnames(catalogue))) {
    stop("Target relations must be scoped by cell type.", call. = FALSE)
  }
  key <- paste(
    as.character(catalogue[[celltype_col]]),
    as.character(catalogue$reaction_id), sep = "\001"
  )
  rows <- split(seq_len(nrow(catalogue)), key)
  answer <- do.call(rbind, lapply(rows, function(index) {
    one <- catalogue[index, , drop = FALSE]
    stages <- sort(unique(as.character(
      one$merged_catalogue_inclusion_stage[
        !is.na(one$merged_catalogue_inclusion_stage) &
          nzchar(one$merged_catalogue_inclusion_stage)
      ]
    )))
    out <- data.frame(
      catalogue_id = as.character(one$catalogue_id[[1L]]),
      cell_type = as.character(one[[celltype_col]][[1L]]),
      reaction_id = as.character(one$reaction_id[[1L]]),
      anchor_core_reaction_ids = paste(
        sort(unique(as.character(one$anchor_core_reaction_id))),
        collapse = ";"
      ),
      expansion_types = paste(
        sort(unique(as.character(one$expansion_type))), collapse = ";"
      ),
      source_annotations = paste(
        sort(unique(as.character(one$source_annotation))), collapse = ";"
      ),
      merged_catalogue_is_core = FALSE,
      merged_catalogue_inclusion_stage = if (length(stages)) {
        paste(stages, collapse = ";")
      } else {
        NA_character_
      },
      score_target = TRUE,
      target_role = "direct_database_crossref_noncore",
      lp_exclusion_reason = NA_character_,
      stringsAsFactors = FALSE
    )
    names(out)[names(out) == "cell_type"] <- celltype_col
    out
  }))
  rownames(answer) <- NULL
  answer
}
'''
replace_function(target_path, '.rc_target_union_aggregate_targets', aggregate_targets)

model_summary = r'''.rc_target_union_model_summary <- function(layer2) {
  if (!inherits(layer2, "regcompass_layer2_step") ||
      !identical(as.character(layer2$model_mode), "meta_module_gem")) {
    stop(
      "`layer2` must be a completed cell-type meta-module GEM run.",
      call. = FALSE
    )
  }
  summary <- layer2$model_cache_summary
  required <- c(
    "cell_type", "medium_scenario", "file", "file_checksum",
    "build_strategy", "completion_stage"
  )
  if (!is.data.frame(summary) || !all(required %in% colnames(summary)) ||
      !nrow(summary)) {
    stop("Layer 2 does not identify cell-type union GEM files.",
         call. = FALSE)
  }
  summary$cell_type <- trimws(as.character(summary$cell_type))
  summary$medium_scenario <- trimws(as.character(summary$medium_scenario))
  summary$file <- as.character(summary$file)
  summary$file_checksum <- as.character(summary$file_checksum)
  valid <- !is.na(summary$cell_type) & nzchar(summary$cell_type) &
    !is.na(summary$medium_scenario) & nzchar(summary$medium_scenario) &
    !is.na(summary$file) & nzchar(summary$file)
  summary <- unique(summary[valid, , drop = FALSE])
  key <- paste(summary$cell_type, summary$medium_scenario, sep = "\001")
  files <- split(summary$file, key)
  ambiguous <- names(files)[vapply(
    files, function(x) length(unique(x)) != 1L, logical(1)
  )]
  if (length(ambiguous)) {
    stop("Each cell type and medium must resolve to one union GEM file.",
         call. = FALSE)
  }
  summary <- summary[!duplicated(key), , drop = FALSE]
  for (i in seq_len(nrow(summary))) {
    if (!identical(summary$build_strategy[[i]],
                   "celltype_medium_union_gem") ||
        !identical(summary$completion_stage[[i]],
                   "celltype_specific_fastcore_after_condition_module_union")) {
      stop("Target remapping requires cell-type-specific final union GEMs.",
           call. = FALSE)
    }
    .rc_read_celltype_union_gem(
      summary$file[[i]], summary$cell_type[[i]],
      summary$medium_scenario[[i]], summary$file_checksum[[i]]
    )
  }
  summary[order(summary$cell_type, summary$medium_scenario), , drop = FALSE]
}
'''
replace_function(target_path, '.rc_target_union_model_summary', model_summary)

cached_ids = r'''.rc_target_union_cached_reaction_ids <- function(layer2) {
  summary <- .rc_target_union_model_summary(layer2)
  by_celltype <- split(seq_len(nrow(summary)), summary$cell_type)
  available <- lapply(by_celltype, function(rows) {
    sets <- lapply(rows, function(i) {
      model <- .rc_read_celltype_union_gem(
        summary$file[[i]], summary$cell_type[[i]],
        summary$medium_scenario[[i]], summary$file_checksum[[i]]
      )
      rc_validate_gem(model)$reactions
    })
    value <- .rc_target_union_normalize_ids(Reduce(intersect, sets))
    if (!length(value)) {
      stop(
        "Cell type `", summary$cell_type[[rows[[1L]]]],
        "` has no reactions shared across its medium-specific union GEMs.",
        call. = FALSE
      )
    }
    value
  })
  attr(available, "model_summary") <- summary
  available
}
'''
replace_function(target_path, '.rc_target_union_cached_reaction_ids', cached_ids)

build_target_cache = r'''.rc_build_target_union_model_cache <- function(
    layer2, target_reactions,
    target_direction = c("both", "forward", "reverse"),
    celltype_col = "cell_type") {
  target_direction <- match.arg(target_direction)
  required <- c(celltype_col, "reaction_id")
  if (!is.data.frame(target_reactions) ||
      !all(required %in% colnames(target_reactions)) || !nrow(target_reactions)) {
    stop("Target reactions must be a non-empty cell-type-scoped table.",
         call. = FALSE)
  }
  summary <- .rc_target_union_model_summary(layer2)
  cache <- list()
  diagnostics <- list()
  fingerprints <- character(nrow(summary))
  for (i in seq_len(nrow(summary))) {
    cell_type <- summary$cell_type[[i]]
    scenario <- summary$medium_scenario[[i]]
    file <- summary$file[[i]]
    checksum <- summary$file_checksum[[i]]
    model <- .rc_read_celltype_union_gem(
      file, cell_type, scenario, checksum
    )
    validated <- rc_validate_gem(model)
    targets <- unique(as.character(target_reactions$reaction_id[
      as.character(target_reactions[[celltype_col]]) == cell_type
    ]))
    if (!length(targets)) next
    fingerprints[[i]] <- .rc_full_gem_cache_fingerprint(model)
    missing <- setdiff(targets, validated$reactions)
    if (length(missing)) {
      stop(
        "Targeted reactions are absent from the `", cell_type,
        "` union GEM for medium `", scenario, "`: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    directions <- rc_prepare_directional_targets(
      model, targets, target_direction
    )
    directions$cell_type <- cell_type
    directions$medium_scenario <- scenario
    diagnostics[[paste(cell_type, scenario, sep = "\001")]] <- directions
    allowed <- directions[
      directions$target_direction %in% c("forward", "reverse"),
      , drop = FALSE
    ]
    for (j in seq_len(nrow(allowed))) {
      reaction <- as.character(allowed$reaction_id[[j]])
      direction <- as.character(allowed$target_direction[[j]])
      key <- paste0(
        "celltype=", utils::URLencode(cell_type, reserved = TRUE),
        "::reaction=", utils::URLencode(reaction, reserved = TRUE),
        "::direction=", direction,
        "::medium=", utils::URLencode(scenario, reserved = TRUE)
      )
      cache[[key]] <- list(
        module_id = "CELLTYPE_MEDIUM_UNION_GEM",
        cell_type = cell_type,
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        file_checksum = checksum,
        build_strategy = "reuse_exact_celltype_medium_union_gem"
      )
    }
  }
  if (!length(cache)) {
    stop("No cell-type targeted reaction direction can be scored.",
         call. = FALSE)
  }
  summary$source_model_fingerprint <- fingerprints
  summary$reused_without_rebuilding <- TRUE
  summary$second_pass_build_strategy <-
    "reuse_exact_celltype_medium_union_gem"
  attr(cache, "summary") <- summary
  attr(cache, "direction_diagnostics") <- .rc_bind_frames_fill(diagnostics)
  attr(cache, "structural_scope") <- "cell_type_x_medium"
  cache
}
'''
replace_function(target_path, '.rc_build_target_union_model_cache', build_target_cache)

score_existing = r'''.rc_score_existing_union_cache <- function(
    layer1, gem, model_cache,
    condition_col, celltype_col,
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    parallel = TRUE, BPPARAM = NULL) {
  summary <- attr(model_cache, "summary")
  if (!is.data.frame(summary) ||
      !all(c("cell_type", "medium_scenario") %in% colnames(summary))) {
    stop("Targeted model cache lacks cell-type provenance.", call. = FALSE)
  }
  media <- unique(as.character(summary$medium_scenario))
  medium_scenarios <- data.frame(
    medium_scenario_id = media,
    exchange_reaction_id = NA_character_,
    lb = NA_real_, ub = NA_real_, available = FALSE,
    .no_constraints = TRUE,
    stringsAsFactors = FALSE
  )
  answer <- .rc_run_celltype_microcompass_engine(
    layer1 = layer1,
    gem = gem,
    target_reactions = NULL,
    medium_scenarios = medium_scenarios,
    reaction_membership = NULL,
    core_reactions = NULL,
    unit = "metacell",
    condition_col = condition_col,
    sample_col = NULL,
    celltype_col = celltype_col,
    model_params = list(),
    omega = omega,
    target_direction = "both",
    parallel = parallel,
    solver = match.arg(solver),
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM,
    model_cache_override = model_cache
  )
  answer$model_mode <- "reused_celltype_medium_union_gem"
  answer$params$structural_model_reused_exactly <- TRUE
  answer$params$fastcore_rerun <- FALSE
  answer$params$model_rebuild <- FALSE
  answer$params$parallel_task <-
    "reused_celltype_union_gem_by_matching_metacell"
  answer$direction_diagnostics <- attr(model_cache, "direction_diagnostics")
  answer$relative_penalty_rank <- answer$score
  answer$score_semantics <- attr(answer$score, "score_semantics") %||%
    "within_target_relative_penalty_rank_not_probability"
  answer$noninformative_target <- attr(answer$score, "noninformative_target")
  answer$primary_output <- "penalty"
  answer$primary_output_semantics <-
    "minimum evidence-discordance penalty; lower means stronger support"
  class(answer) <- c("regcompass_expanded_layer2_result", "list")
  answer
}
'''
replace_function(target_path, '.rc_score_existing_union_cache', score_existing)

# Public targeted remapping validates scoped source rows and passes scoped targets.
replace_once(
    target_path,
    '  if (!setequal(\n    as.character(layer2$source_core_reactions$reaction_id),\n    as.character(catalogue$merged_core_reactions$reaction_id)\n  )) {\n    stop("Layer 2 was not generated from the supplied merged core reactions.",\n         call. = FALSE)\n  }\n',
    '  source_key <- paste(\n'
    '    as.character(layer2$source_core_reactions[[workflow$celltype_col]]),\n'
    '    as.character(layer2$source_core_reactions$reaction_id), sep = "\\001"\n'
    '  )\n'
    '  catalogue_key <- paste(\n'
    '    as.character(catalogue$merged_core_reactions[[workflow$celltype_col]]),\n'
    '    as.character(catalogue$merged_core_reactions$reaction_id), sep = "\\001"\n'
    '  )\n'
    '  if (!setequal(source_key, catalogue_key)) {\n'
    '    stop("Layer 2 was not generated from the supplied cell-type cores.",\n'
    '         call. = FALSE)\n'
    '  }\n'
)
replace_once(
    target_path,
    '    cached_reaction_ids = cached_reaction_ids\n  )\n',
    '    cached_reaction_ids = cached_reaction_ids,\n'
    '    celltype_col = workflow$celltype_col\n'
    '  )\n'
)
replace_once(
    target_path,
    '    target_reactions = definition$params$score_targets,\n    target_direction = target_direction\n  )\n',
    '    target_reactions = definition$params$score_targets,\n'
    '    target_direction = target_direction,\n'
    '    celltype_col = workflow$celltype_col\n'
    '  )\n'
)
replace_once(
    target_path,
    '  scored$params$n_expanded_score_targets <-\n    length(definition$params$score_targets)\n',
    '  scored$params$n_expanded_score_targets <-\n'
    '    nrow(definition$params$score_targets)\n'
)

# Update architecture tests.
(TESTS / 'test_union_gem_architecture.R').write_text(r'''test_that("meta-modules merge conditions within cell type only", {
  biological <- data.frame(
    group_id = c("C1|T", "C2|T", "C1|B"),
    condition = c("C1", "C2", "C1"),
    cell_type = c("T", "T", "B"),
    module_id = c("T1", "T2", "B1"),
    reaction_id = c("RT1", "RT2", "RB1"),
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  condition_modules <- list(
    condition_fit_status = biological[, c("group_id", "condition", "cell_type")],
    tf_peak_gene_condition_all = data.frame(),
    tf_peak_gene_condition = data.frame(),
    supported_metabolic_genes = data.frame(),
    core_gene_reaction = biological,
    reaction_membership = biological,
    meta_module_summary = data.frame()
  )
  merged <- .rc_merge_meta_modules_by_cell_type(
    condition_modules, "cell_type", "condition"
  )
  expect_setequal(names(merged$cell_type_catalogues), c("T", "B"))
  expect_setequal(
    merged$cell_type_catalogues$T$merged_core_reactions$reaction_id,
    c("RT1", "RT2")
  )
  expect_identical(
    merged$cell_type_catalogues$B$merged_core_reactions$reaction_id,
    "RB1"
  )
  expect_true(all(
    merged$merged_core_reactions$cell_type ==
      c("B", "T", "T")[match(
        merged$merged_core_reactions$reaction_id,
        c("RB1", "RT1", "RT2")
      )]
  ))
  expect_identical(merged$merge_scope, "cell_type")
  expect_false(merged$cross_celltype_merge)
  expect_false(merged$is_gem)
  expect_false(merged$fastcore_applied)
})

test_that("Stage 3 contains no FASTCORE execution path", {
  construction <- paste(
    deparse(body(.rc_build_condition_meta_modules)), collapse = "\n"
  )
  stage <- paste(deparse(body(rc_regcompass_step_meta_modules)), collapse = "\n")
  expect_false(grepl(".rc_complete_celltype_medium_union_gem", construction,
                     fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", construction, fixed = TRUE))
  expect_false(grepl(".rc_complete_celltype_medium_union_gem", stage,
                     fixed = TRUE))
  expect_false(grepl(".rc_fastcore_", stage, fixed = TRUE))
  expect_true(grepl("none_at_meta_module_stage", construction, fixed = TRUE))
  expect_true(grepl("merge_creates_gem = FALSE", stage, fixed = TRUE))
})

test_that("union GEM and FASTCORE scopes are cell type by medium", {
  cache_body <- paste(
    deparse(body(.rc_build_celltype_medium_union_gem_cache)), collapse = "\n"
  )
  completion_body <- paste(
    deparse(body(.rc_complete_celltype_medium_union_gem)), collapse = "\n"
  )
  engine_body <- paste(
    deparse(body(.rc_run_celltype_microcompass_engine)), collapse = "\n"
  )
  expect_true(grepl("for (cell_type in scoped$cell_types)", cache_body,
                    fixed = TRUE))
  expect_true(grepl("celltype_specific_fastcore", cache_body, fixed = TRUE))
  expect_true(grepl("cell_type = cell_type", completion_body, fixed = TRUE))
  expect_true(grepl("celltype_fastcore_support", completion_body, fixed = TRUE))
  expect_true(grepl("unit_celltype == cell_type", engine_body, fixed = TRUE))
  expect_true(grepl("shared_across_cell_types = FALSE", engine_body,
                    fixed = TRUE))
})

test_that("microCOMPASS row IDs retain cell-type structural scope", {
  parsed <- rc_parse_microcompass_row_id(
    paste0(
      "celltype=Tumor%20cell::reaction=R1::direction=forward::medium=base"
    )
  )
  expect_identical(parsed$cell_type, "Tumor cell")
  expect_identical(parsed$reaction_id, "R1")
})
''')

# Update FASTCORE tests and all static references.
fast_test = TESTS / 'test_fastcore_architecture.R'
text = fast_test.read_text()
text = text.replace('.rc_complete_medium_union_gem(',
                    '.rc_complete_celltype_medium_union_gem(')
# Every test call uses a named gem argument, so insert a deterministic scope.
text = text.replace(
    '    gem = gem,\n    reaction_membership =',
    '    gem = gem,\n    cell_type = "test_cell_type",\n    reaction_membership ='
)
text = text.replace('global add-only FASTCORE', 'cell-type add-only FASTCORE')
text = text.replace('global_fastcore_support', 'celltype_fastcore_support')
text = text.replace(
    '"single_global_fastcore_after_meta_module_merge"',
    '"celltype_specific_fastcore_after_condition_module_union"'
)
text = text.replace(
    'sample=S1::module=S1%3A%3AM1::reaction=R1::direction=forward::medium=base',
    'sample=S1::module=S1%3A%3AM1::celltype=T::reaction=R1::direction=forward::medium=base'
)
text = text.replace(
    '  expect_equal(parsed$module_id, "S1::M1")\n',
    '  expect_equal(parsed$module_id, "S1::M1")\n'
    '  expect_equal(parsed$cell_type, "T")\n'
)
fast_test.write_text(text)

linux_test = TESTS / 'test-linux-parallel-workflow.R'
linux_text = linux_test.read_text().replace(
    '.rc_complete_medium_union_gem(',
    '.rc_complete_celltype_medium_union_gem('
)
linux_test.write_text(linux_text)

pando_test = TESTS / 'test_pando_meta_module.R'
ptext = pando_test.read_text()
ptext = ptext.replace(
    '    group_id = "S1",\n    sample_id = "S1",\n',
    '    group_id = "S1",\n    sample_id = "S1",\n    condition = "C1",\n    cell_type = "T",\n',
    1
)
ptext = ptext.replace(
    '    sample_id = "S1",\n    module_id = "S1::SUPPORTED_METABOLIC_GENES",\n    reaction_id = c("R1", "R2"),\n',
    '    sample_id = "S1",\n    condition = "C1",\n    cell_type = "T",\n    module_id = "S1::SUPPORTED_METABOLIC_GENES",\n    reaction_id = c("R1", "R2"),\n',
    1
)
ptext = ptext.replace(
    '    condition_fit_status = data.frame(group_id = "S1"),',
    '    condition_fit_status = data.frame(\n'
    '      group_id = "S1", condition = "C1", cell_type = "T"\n'
    '    ),'
)
ptext = ptext.replace(
    '  merged <- .rc_merge_meta_module_catalogue(condition_modules)',
    '  merged <- .rc_merge_meta_modules_by_cell_type(\n'
    '    condition_modules, "cell_type", "condition"\n'
    '  )'
)
ptext = ptext.replace(
    '  union_gem <- .rc_complete_medium_union_gem(\n    gem = gem,',
    '  union_gem <- .rc_complete_celltype_medium_union_gem(\n'
    '    gem = gem,\n    cell_type = "T",'
)
ptext = ptext.replace('global_fastcore_support', 'celltype_fastcore_support')
pando_test.write_text(ptext)

# Documentation: explicit structural boundary and terminology.
appendices = {
    ROOT / 'docs' / 'tutorial-02-stepwise-audit.md': r'''

## Cell-type boundary for meta-modules and structural models

Stage 3 retains condition-specific evidence, but biological meta-modules are
unioned **only within the same cell type**. A reaction supported in one cell
type is never inserted into another cell type's catalogue.

Stage 5 therefore constructs one structural model for every
`cell_type × medium_scenario` combination. Conditions and metacells from the
same cell type reuse that model. Different cell types have separate reaction
memberships, separate union GEM files, separate FASTCORE completion runs,
separate model checksums, and separate directional `vmax` caches.
''',
    ROOT / 'docs' / 'tutorial-03-mathematical-model.md': r'''

## Cell-type-specific union GEM

Let `c` denote cell type, `d` condition, and `m` medium. Condition-specific
biological reaction sets are unioned only within cell type,
`B_c = union_d B_{c,d}`. RegCompass then runs FASTCORE independently for every
`(c,m)` pair and obtains `G_{c,m}`. There is no operation of the form
`union_c G_{c,m}`. Consequently, conditions are structurally comparable within
a cell type, while biologically distinct cell types are not forced onto an
artificial cross-cell-type reaction universe.
''',
    ROOT / 'docs' / 'tutorial-04-targeted-reaction-remapping.md': r'''

## Cell-type-scoped cache reuse

Targeted remapping preserves the original structural boundary. A target linked
to a core reaction from cell type `c` is evaluated only in cached
`cell_type = c` union GEMs. The function intersects candidate reactions across
media within that cell type; it never intersects or merges cached GEMs from
different cell types and never reruns FASTCORE.
''',
    ROOT / 'docs' / 'tutorial-05-condition-differential-analysis.md': r'''

## Valid comparison scope

Condition contrasts are performed for the same reaction, direction, medium,
and cell type. Each row is evaluated only in metacells whose cell type matches
the row's union-GEM scope. Cross-cell-type penalty or `vmax` comparisons are
not treated as condition contrasts.
'''
}
for path, addition in appendices.items():
    text = path.read_text()
    marker = addition.strip().splitlines()[0]
    if marker not in text:
        path.write_text(text.rstrip() + '\n' + addition)

news = ROOT / 'NEWS.md'
if news.exists():
    text = news.read_text()
    entry = (
        '- Union meta-modules, union GEMs, FASTCORE completion, model caches, '
        'and directional vmax reuse are now partitioned by cell type; only '
        'conditions within the same cell type share a structural model.\n'
    )
    if entry not in text:
        news.write_text(entry + text)

# Eliminate misleading global terminology in current tutorials and source docs.
for path in list((ROOT / 'docs').glob('*.md')) + [ROOT / 'DESCRIPTION']:
    text = path.read_text()
    text = text.replace(
        'single global FASTCORE after meta-module merge',
        'independent FASTCORE within each cell-type union GEM'
    )
    text = text.replace(
        'one final union GEM per medium shared across all units',
        'one union GEM per cell type and medium shared only within that cell type'
    )
    text = text.replace(
        'one medium-specific union GEM',
        'one cell-type-specific union GEM per medium'
    )
    path.write_text(text)

# Static invariants after all transformations.
all_r = '\n'.join(path.read_text() for path in R.glob('*.R'))
for retired in (
    '.rc_merge_meta_module_catalogue',
    '.rc_build_medium_specific_union_gem_cache',
    '.rc_complete_medium_union_gem',
    '.rc_read_cached_union_gem',
    'one_final_union_gem_per_medium_shared_across_all_units',
    'single_global_fastcore_after_meta_module_merge'
):
    if retired in all_r:
        raise RuntimeError(f'Retired cross-cell-type implementation remains: {retired}')
