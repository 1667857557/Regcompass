# Score direct database-linked non-core reactions in final union GEMs.

.rc_target_union_normalize_ids <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

.rc_target_union_core_rows <- function(
    gem, available_core_reactions,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct")) {
  gene_match <- match.arg(gene_match)
  validated <- rc_validate_gem(gem)
  original_core <- intersect(
    .rc_target_union_normalize_ids(available_core_reactions),
    validated$reactions
  )
  if (!length(original_core)) {
    stop("The original Layer 2 run contains no valid core reactions.",
         call. = FALSE)
  }

  requested_reactions <- .rc_target_union_normalize_ids(core_reaction_ids)
  requested_genes <- toupper(.rc_target_union_normalize_ids(core_genes))
  if (!length(requested_reactions) && !length(requested_genes)) {
    stop(
      "Supply at least one `core_reaction_ids` or `core_genes` value.",
      call. = FALSE
    )
  }

  missing_reactions <- setdiff(requested_reactions, validated$reactions)
  if (length(missing_reactions)) {
    stop(
      "Selected reactions are absent from the GEM: ",
      paste(utils::head(missing_reactions, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  gpr <- gem$gpr_table
  gene_reactions <- character()
  gene_source <- list()
  if (length(requested_genes)) {
    required <- c("reaction_id", "and_group_id", "gene")
    if (!is.data.frame(gpr) || !all(required %in% colnames(gpr))) {
      stop(
        paste(
          "Gene-selected cores require a GEM `gpr_table` containing",
          "reaction_id, and_group_id and gene."
        ),
        call. = FALSE
      )
    }
    gpr <- unique(gpr[, required, drop = FALSE])
    gpr$reaction_id <- trimws(as.character(gpr$reaction_id))
    gpr$and_group_id <- as.character(gpr$and_group_id)
    gpr$gene <- toupper(trimws(as.character(gpr$gene)))
    gpr <- gpr[
      gpr$reaction_id %in% validated$reactions &
        !is.na(gpr$gene) & nzchar(gpr$gene),
      , drop = FALSE
    ]
    absent_genes <- setdiff(requested_genes, unique(gpr$gene))
    if (length(absent_genes)) {
      stop(
        "Selected core genes do not map to GEM GPR rules: ",
        paste(utils::head(absent_genes, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    if (identical(gene_match, "any_direct")) {
      mapped <- unique(gpr$reaction_id[gpr$gene %in% requested_genes])
    } else {
      group_key <- paste(gpr$reaction_id, gpr$and_group_id, sep = "\001")
      groups <- split(seq_len(nrow(gpr)), group_key)
      complete <- vapply(groups, function(rows) {
        all(unique(gpr$gene[rows]) %in% requested_genes)
      }, logical(1))
      mapped <- unique(vapply(groups[complete], function(rows) {
        gpr$reaction_id[rows[[1L]]]
      }, character(1)))
    }
    gene_reactions <- intersect(mapped, original_core)
    if (!length(gene_reactions) && !length(requested_reactions)) {
      stop(
        "The selected genes do not resolve to original Layer 2 core targets.",
        call. = FALSE
      )
    }
    gene_source <- lapply(gene_reactions, function(reaction) {
      sort(unique(gpr$gene[
        gpr$reaction_id == reaction & gpr$gene %in% requested_genes
      ]))
    })
    names(gene_source) <- gene_reactions
  }

  reactions <- union(requested_reactions, gene_reactions)
  if (!length(reactions)) {
    stop("No selected anchor is available in this cell type.",
         call. = FALSE)
  }
  selection_source <- vapply(reactions, function(reaction) {
    by_id <- reaction %in% requested_reactions
    by_gene <- reaction %in% gene_reactions
    if (by_id && by_gene) {
      "reaction_id+core_gene"
    } else if (by_id) {
      "reaction_id"
    } else if (identical(gene_match, "complete_gpr")) {
      "core_gene_complete_gpr"
    } else {
      "core_gene_any_direct"
    }
  }, character(1))
  mapped_genes <- vapply(reactions, function(reaction) {
    genes <- gene_source[[reaction]]
    if (is.null(genes) || !length(genes)) NA_character_ else
      paste(genes, collapse = ";")
  }, character(1))
  is_core <- reactions %in% original_core

  data.frame(
    catalogue_id = "MERGED_META_MODULES",
    gene = mapped_genes,
    reaction_id = reactions,
    is_core = is_core,
    anchor_role = ifelse(is_core, "original_layer2_core", "gem_noncore"),
    selection_source = selection_source,
    stringsAsFactors = FALSE
  )
}

.rc_build_target_union_definition <- function(
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

.rc_target_union_direct_crossref_relations <- function(
    gem, selected_core_reactions) {
  maps <- rc_reaction_crossref_maps(gem)
  specifications <- list(
    list(
      map = .rc_clean_meta_module_map(maps$kegg, "kegg_id"),
      id_col = "kegg_id",
      expansion_type = "shared_kegg_reaction",
      prefix = "KEGG:"
    ),
    list(
      map = .rc_clean_meta_module_map(maps$reactome, "reactome_id"),
      id_col = "reactome_id",
      expansion_type = "shared_reactome_reaction",
      prefix = "REACTOME:"
    ),
    list(
      map = .rc_clean_meta_module_map(
        maps$rhea_master, "rhea_master_id"
      ),
      id_col = "rhea_master_id",
      expansion_type = "shared_master_rhea_reaction",
      prefix = "RHEA_MASTER:"
    )
  )
  anchors <- .rc_target_union_normalize_ids(
    selected_core_reactions$reaction_id
  )
  output <- list()
  output_index <- 0L
  for (anchor in anchors) {
    for (specification in specifications) {
      map <- specification$map
      id_col <- specification$id_col
      if (!is.data.frame(map) || !nrow(map)) next
      anchor_ids <- unique(as.character(map[[id_col]][
        map$reaction_id == anchor
      ]))
      anchor_ids <- anchor_ids[
        !is.na(anchor_ids) & nzchar(trimws(anchor_ids))
      ]
      if (!length(anchor_ids)) next
      reactions <- unique(as.character(map$reaction_id[
        map[[id_col]] %in% anchor_ids
      ]))
      reactions <- setdiff(reactions, anchor)
      for (reaction in reactions) {
        shared_ids <- intersect(
          anchor_ids,
          unique(as.character(map[[id_col]][map$reaction_id == reaction]))
        )
        shared_ids <- sort(shared_ids[
          !is.na(shared_ids) & nzchar(trimws(shared_ids))
        ])
        if (!length(shared_ids)) next
        output_index <- output_index + 1L
        output[[output_index]] <- data.frame(
          anchor_core_reaction_id = anchor,
          reaction_id = reaction,
          expansion_type = specification$expansion_type,
          source_annotation = paste0(
            specification$prefix, paste(shared_ids, collapse = ";")
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(output)) {
    return(data.frame(
      anchor_core_reaction_id = character(),
      reaction_id = character(),
      expansion_type = character(),
      source_annotation = character(),
      stringsAsFactors = FALSE
    ))
  }
  answer <- unique(do.call(rbind, output))
  answer <- answer[order(
    answer$anchor_core_reaction_id,
    answer$reaction_id,
    answer$expansion_type,
    answer$source_annotation
  ), , drop = FALSE]
  rownames(answer) <- NULL
  answer
}

.rc_target_union_aggregate_targets <- function(
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

.rc_target_union_model_summary <- function(layer2) {
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

.rc_target_union_cached_reaction_ids <- function(layer2) {
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

.rc_build_target_union_model_cache <- function(
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

.rc_score_existing_union_cache <- function(
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

#' Score directly database-linked non-core reactions in final union GEMs
#'
#' Selected original core reactions are mapping anchors only. The function
#' scores directly KEGG-, Reactome-, or master-Rhea-linked non-core reactions
#' by reusing the exact final medium-specific union GEM files created by Stage
#' 5. It does not rebuild a GEM and does not rerun FASTCORE.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer2 Output from [rc_regcompass_step_layer2()] with
#'   `model_mode = "meta_module_gem"`.
#' @param gem The same GEM used for the original run.
#' @param outdir Output directory.
#' @param core_reaction_ids GEM reaction IDs used as direct mapping anchors.
#'   The historical argument name is retained for compatibility; an ID may be
#'   an original Layer 2 core or any other valid reaction in `gem`.
#' @param core_genes Genes used to resolve original core anchors through GPRs.
#'   Gene selection intentionally remains restricted to original cores.
#' @param gene_match Require a complete GPR group or allow any direct gene match.
#' @param layer2_args Optional `omega`, `target_direction`, `solver`, and
#'   `flux_threshold` overrides. Scoring LPs have no time-limit control.
#' @param parallel Whether to parallelize model-by-metacell tasks.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param progress Whether to display stage progress.
#' @return A `regcompass_target_union_step` containing selected anchors, direct
#'   database-linked non-core targets, final union-GEM provenance, and LP scores.
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("target_union", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gene_match <- match.arg(gene_match)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  workflow <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = workflow, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = workflow, gem = gem,
    required_mode = "meta_module_gem", argument = "layer2"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c("omega", "target_direction", "solver", "flux_threshold")
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop(
      "Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
      ". Scoring `time_limit` has been removed; only ",
      "`layer2_args$model_params$completion_time_limit` in the original ",
      "union-GEM construction stage is supported.",
      call. = FALSE
    )
  }
  catalogue <- meta_modules$merged_modules
  if (!is.list(catalogue) ||
      !is.data.frame(catalogue$merged_core_reactions) ||
      !is.data.frame(catalogue$merged_reaction_membership)) {
    stop("The merged biological meta-module catalogue is unavailable.",
         call. = FALSE)
  }
  source_key <- paste(
    as.character(layer2$source_core_reactions[[workflow$celltype_col]]),
    as.character(layer2$source_core_reactions$reaction_id), sep = "\001"
  )
  catalogue_key <- paste(
    as.character(catalogue$merged_core_reactions[[workflow$celltype_col]]),
    as.character(catalogue$merged_core_reactions$reaction_id), sep = "\001"
  )
  if (!setequal(source_key, catalogue_key)) {
    stop("Layer 2 was not generated from the supplied cell-type cores.",
         call. = FALSE)
  }
  cached_reaction_ids <- .rc_target_union_cached_reaction_ids(layer2)
  definition <- .rc_build_target_union_definition(
    gem = gem,
    merged_core_reactions = catalogue$merged_core_reactions,
    merged_reaction_membership = catalogue$merged_reaction_membership,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    cached_reaction_ids = cached_reaction_ids,
    celltype_col = workflow$celltype_col
  )
  target_direction <- match.arg(
    as.character(layer2_args$target_direction %||%
                   layer2$params$target_direction %||% "both"),
    c("both", "forward", "reverse")
  )
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  omega <- layer2_args$omega %||% layer2$params$omega %||% 0.95
  flux_threshold <- layer2_args$flux_threshold %||% 1e-8
  model_cache <- .rc_build_target_union_model_cache(
    layer2 = layer2,
    target_reactions = definition$params$score_targets,
    target_direction = target_direction,
    celltype_col = workflow$celltype_col
  )
  scored <- .rc_score_existing_union_cache(
    layer1 = layer1,
    gem = gem,
    model_cache = model_cache,
    condition_col = workflow$condition_col,
    celltype_col = workflow$celltype_col,
    omega = omega,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  scored$workflow_params <- workflow
  scored$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  scored$params$target_direction <- target_direction
  scored$params$target_scope <-
    "direct_kegg_reactome_master_rhea_noncore_only"
  scored$params$n_selected_core <- nrow(definition$selected_core_reactions)
  scored$params$n_merged_core_reactions_not_rescored <-
    length(definition$params$merged_core_reactions_not_rescored)
  scored$params$n_cached_union_unavailable_reactions <-
    length(definition$params$cached_union_unavailable_reactions)
  scored$params$n_expanded_score_targets <-
    nrow(definition$params$score_targets)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_mm_write_tsv_gz(
    definition$selected_core_reactions,
    file.path(outdir, "selected_core_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$expanded_reaction_catalog,
    file.path(outdir, "expanded_reaction_catalog.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$expanded_scoring_targets,
    file.path(outdir, "expanded_scoring_targets.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$merged_catalogue_membership,
    file.path(outdir, "merged_meta_module_catalogue_membership.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$summary,
    file.path(outdir, "target_union_summary.tsv.gz")
  )
  rc_export_microcompass(scored, file.path(outdir, "scores"))
  answer <- c(definition, list(microcompass = scored))
  answer$workflow_params <- workflow
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(answer) <- c("regcompass_target_union_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
