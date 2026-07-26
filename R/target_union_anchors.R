# Extend target-union remapping to any valid GEM reaction anchor while
# preserving the existing public argument names and core-only output fields.

.rc_target_union_definition_core_only <- .rc_build_target_union_definition
.rc_target_union_step_core_only <- rc_regcompass_step_target_union

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
          "Gene-selected anchors require a GEM `gpr_table` containing",
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
        "Selected genes do not map to GEM GPR rules: ",
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
    gene_reactions <- intersect(mapped, validated$reactions)
    if (!length(gene_reactions)) {
      stop(
        "The selected genes do not resolve to GEM reaction anchors under `gene_match`.",
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
      "reaction_id+gene"
    } else if (by_id) {
      "reaction_id"
    } else if (identical(gene_match, "complete_gpr")) {
      "gene_complete_gpr"
    } else {
      "gene_any_direct"
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
  answer <- tryCatch(
    .rc_target_union_definition_core_only(
      gem = gem,
      merged_core_reactions = merged_core_reactions,
      merged_reaction_membership = merged_reaction_membership,
      core_reaction_ids = core_reaction_ids,
      core_genes = core_genes,
      gene_match = gene_match,
      cached_reaction_ids = cached_reaction_ids
    ),
    error = function(error) {
      message <- conditionMessage(error)
      message <- sub(
        "The selected core reactions have no directly linked",
        "The selected anchor reactions have no directly linked",
        message,
        fixed = TRUE
      )
      stop(message, call. = FALSE)
    }
  )

  anchors <- answer$selected_core_reactions
  if (!"is_core" %in% colnames(anchors)) {
    anchors$is_core <- anchors$reaction_id %in%
      as.character(merged_core_reactions$reaction_id)
  }
  if (!"anchor_role" %in% colnames(anchors)) {
    anchors$anchor_role <- ifelse(
      anchors$is_core, "original_layer2_core", "gem_noncore"
    )
  }
  core_anchors <- anchors[anchors$is_core %in% TRUE, , drop = FALSE]
  noncore_anchors <- anchors[!anchors$is_core %in% TRUE, , drop = FALSE]

  catalogue <- answer$expanded_reaction_catalog
  catalogue$anchor_reaction_id <- catalogue$anchor_core_reaction_id
  anchor_is_core <- stats::setNames(
    as.logical(anchors$is_core), as.character(anchors$reaction_id)
  )
  catalogue$anchor_is_original_core <- unname(
    anchor_is_core[as.character(catalogue$anchor_reaction_id)]
  )
  catalogue$anchor_is_original_core[
    is.na(catalogue$anchor_is_original_core)
  ] <- FALSE

  targets <- answer$expanded_scoring_targets
  if ("anchor_core_reaction_ids" %in% colnames(targets)) {
    targets$anchor_reaction_ids <- targets$anchor_core_reaction_ids
  }

  answer$selected_anchor_reactions <- anchors
  answer$selected_core_reactions <- core_anchors
  answer$selected_noncore_reactions <- noncore_anchors
  answer$expanded_reaction_catalog <- catalogue
  answer$expanded_scoring_targets <- targets
  answer$summary$n_selected_anchors <- nrow(anchors)
  answer$summary$n_selected_core <- nrow(core_anchors)
  answer$summary$n_selected_noncore_anchors <- nrow(noncore_anchors)
  answer$summary$expansion_policy <-
    "direct_from_selected_anchors_via_kegg_reactome_master_rhea_only"
  answer$params$selected_anchor_reactions <- unique(
    as.character(anchors$reaction_id)
  )
  answer$params$selected_core_reactions <- unique(
    as.character(core_anchors$reaction_id)
  )
  answer$params$selected_noncore_reactions <- unique(
    as.character(noncore_anchors$reaction_id)
  )
  answer$params$expansion_policy <-
    "direct_from_selected_anchors_via_kegg_reactome_master_rhea_only"
  answer
}

#' Score database-linked non-core reactions from GEM reaction anchors
#'
#' `core_reaction_ids` and `core_genes` are retained as public argument names for
#' compatibility. They may resolve to either original Layer 2 core reactions or
#' other reactions in the supplied GEM. Anchors are used only for direct KEGG,
#' Reactome, or master-Rhea remapping. A remapped target is scored only when it
#' is non-core and present in every required cached Stage 5 union GEM.
#'
#' @inheritParams rc_regcompass_step_target_union
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  answer <- .rc_target_union_step_core_only(
    layer1 = layer1,
    meta_modules = meta_modules,
    layer2 = layer2,
    gem = gem,
    outdir = outdir,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    layer2_args = layer2_args,
    parallel = parallel,
    BPPARAM = BPPARAM,
    progress = progress
  )

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_mm_write_tsv_gz(
    answer$selected_anchor_reactions,
    file.path(outdir, "selected_anchor_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    answer$selected_noncore_reactions,
    file.path(outdir, "selected_noncore_reactions.tsv.gz")
  )
  answer$microcompass$params$n_selected_anchors <-
    nrow(answer$selected_anchor_reactions)
  answer$microcompass$params$n_selected_core <-
    nrow(answer$selected_core_reactions)
  answer$microcompass$params$n_selected_noncore_anchors <-
    nrow(answer$selected_noncore_reactions)
  answer$microcompass$params$target_scope <-
    "direct_kegg_reactome_master_rhea_noncore_from_any_gem_anchor"
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
