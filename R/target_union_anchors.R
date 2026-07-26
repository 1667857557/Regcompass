# Extend target-union remapping to any valid GEM reaction-ID anchor.
# Gene selectors retain the original core-only behavior.

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
    if (!length(gene_reactions)) {
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
    cached_reaction_ids = NULL) {
  gene_match <- match.arg(gene_match)
  if (!is.data.frame(merged_core_reactions) ||
      !"reaction_id" %in% colnames(merged_core_reactions)) {
    stop("`merged_core_reactions` must contain reaction_id.",
         call. = FALSE)
  }
  if (!is.data.frame(merged_reaction_membership) ||
      !"reaction_id" %in% colnames(merged_reaction_membership)) {
    stop("`merged_reaction_membership` must contain reaction_id.",
         call. = FALSE)
  }

  selected_anchors <- .rc_target_union_core_rows(
    gem = gem,
    available_core_reactions = merged_core_reactions$reaction_id,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
  selected_core <- selected_anchors[
    selected_anchors$is_core %in% TRUE, , drop = FALSE
  ]
  selected_noncore <- selected_anchors[
    !selected_anchors$is_core %in% TRUE, , drop = FALSE
  ]

  catalogue <- .rc_target_union_direct_crossref_relations(
    gem, selected_anchors
  )
  if (!nrow(catalogue)) {
    stop(
      paste(
        "The selected anchor reactions have no directly linked KEGG, Reactome,",
        "or master-Rhea reactions."
      ),
      call. = FALSE
    )
  }

  catalogue$anchor_reaction_id <- catalogue$anchor_core_reaction_id
  anchor_is_core <- stats::setNames(
    as.logical(selected_anchors$is_core),
    as.character(selected_anchors$reaction_id)
  )
  catalogue$anchor_is_original_core <- unname(
    anchor_is_core[as.character(catalogue$anchor_reaction_id)]
  )
  catalogue$anchor_is_original_core[
    is.na(catalogue$anchor_is_original_core)
  ] <- FALSE

  membership_ids <- .rc_target_union_normalize_ids(
    merged_reaction_membership$reaction_id
  )
  available_ids <- .rc_target_union_normalize_ids(cached_reaction_ids)
  if (!length(available_ids)) {
    stop("No reusable reactions were found in the final union GEM cache.",
         call. = FALSE)
  }

  catalogue_match <- match(
    as.character(catalogue$reaction_id),
    as.character(merged_reaction_membership$reaction_id)
  )
  catalogue$available_in_all_cached_union_gems <-
    catalogue$reaction_id %in% available_ids
  catalogue$present_in_merged_catalogue <- !is.na(catalogue_match)
  catalogue$merged_catalogue_is_core <- catalogue$reaction_id %in%
    as.character(merged_core_reactions$reaction_id)
  catalogue$merged_catalogue_inclusion_stage <- if (
    "inclusion_stage" %in% colnames(merged_reaction_membership)
  ) {
    as.character(merged_reaction_membership$inclusion_stage[catalogue_match])
  } else {
    rep(NA_character_, nrow(catalogue))
  }
  support_only <- catalogue$available_in_all_cached_union_gems &
    !catalogue$present_in_merged_catalogue
  catalogue$merged_catalogue_inclusion_stage[support_only] <-
    "global_fastcore_support_not_in_merged_catalogue"
  catalogue$score_target <- !catalogue$merged_catalogue_is_core &
    catalogue$available_in_all_cached_union_gems
  catalogue$target_role <- ifelse(
    catalogue$merged_catalogue_is_core,
    "merged_core_not_rescored",
    ifelse(
      catalogue$available_in_all_cached_union_gems,
      "direct_database_crossref_noncore",
      "direct_database_crossref_absent_from_cached_union_gem"
    )
  )
  catalogue$lp_exclusion_reason <- ifelse(
    catalogue$merged_catalogue_is_core,
    "already_scored_in_original_layer2",
    ifelse(
      catalogue$available_in_all_cached_union_gems,
      NA_character_,
      "absent_from_one_or_more_cached_union_gems"
    )
  )

  target_relations <- catalogue[catalogue$score_target, , drop = FALSE]
  if (!nrow(target_relations)) {
    stop(
      paste(
        "No directly linked non-core reaction is present in every final",
        "medium-specific union GEM. Original core reactions are not recomputed."
      ),
      call. = FALSE
    )
  }
  targets <- .rc_target_union_aggregate_targets(target_relations)
  targets$anchor_reaction_ids <- targets$anchor_core_reaction_ids

  rownames(selected_anchors) <- NULL
  rownames(selected_core) <- NULL
  rownames(selected_noncore) <- NULL
  rownames(catalogue) <- NULL

  expansion_policy <- if (nrow(selected_noncore)) {
    "direct_from_selected_anchors_via_kegg_reactome_master_rhea_only"
  } else {
    "direct_from_selected_core_via_kegg_reactome_master_rhea_only"
  }
  summary <- data.frame(
    n_selected_anchors = nrow(selected_anchors),
    n_selected_core = nrow(selected_core),
    n_selected_noncore_anchors = nrow(selected_noncore),
    n_direct_crossref_relations = nrow(catalogue),
    n_direct_crossref_reactions = length(unique(catalogue$reaction_id)),
    n_merged_core_reactions_not_rescored = length(unique(
      catalogue$reaction_id[catalogue$merged_catalogue_is_core]
    )),
    n_cached_union_unavailable_reactions = length(unique(
      catalogue$reaction_id[!catalogue$available_in_all_cached_union_gems]
    )),
    n_expanded_score_targets = nrow(targets),
    n_merged_catalogue_reactions = length(membership_ids),
    n_reactions_shared_by_cached_union_gems = length(available_ids),
    gene_match = gene_match,
    expansion_policy = expansion_policy,
    scoring_policy =
      "direct_noncore_reactions_present_in_all_final_union_gems",
    model_policy = "reuse_exact_final_medium_specific_union_gems",
    stringsAsFactors = FALSE
  )

  list(
    selected_anchor_reactions = selected_anchors,
    selected_core_reactions = selected_core,
    selected_noncore_reactions = selected_noncore,
    expanded_reaction_catalog = catalogue,
    expanded_scoring_targets = targets,
    merged_catalogue_membership = merged_reaction_membership,
    summary = summary,
    params = list(
      gene_match = gene_match,
      selected_anchor_reactions = unique(as.character(
        selected_anchors$reaction_id
      )),
      selected_core_reactions = unique(as.character(
        selected_core$reaction_id
      )),
      selected_noncore_reactions = unique(as.character(
        selected_noncore$reaction_id
      )),
      merged_core_reactions_not_rescored = unique(as.character(
        catalogue$reaction_id[catalogue$merged_catalogue_is_core]
      )),
      cached_union_unavailable_reactions = unique(as.character(
        catalogue$reaction_id[!catalogue$available_in_all_cached_union_gems]
      )),
      score_targets = unique(as.character(targets$reaction_id)),
      expansion_policy = expansion_policy
    )
  )
}
