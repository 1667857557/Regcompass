# User-selected core reaction union-GEM construction and scoring.

.rc_target_union_normalize_ids <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

.rc_target_union_core_rows <- function(
    gem, core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct")) {
  gene_match <- match.arg(gene_match)
  validated <- rc_validate_gem(gem)
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
      "Selected core reactions are absent from the GEM: ",
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
        "Gene-selected cores require a GEM `gpr_table` containing reaction_id, and_group_id and gene.",
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
      gene_reactions <- unique(gpr$reaction_id[gpr$gene %in% requested_genes])
    } else {
      group_key <- paste(gpr$reaction_id, gpr$and_group_id, sep = "\001")
      groups <- split(seq_len(nrow(gpr)), group_key)
      complete <- vapply(groups, function(rows) {
        all(unique(gpr$gene[rows]) %in% requested_genes)
      }, logical(1))
      complete_groups <- groups[complete]
      gene_reactions <- unique(vapply(complete_groups, function(rows) {
        gpr$reaction_id[rows[[1L]]]
      }, character(1)))
      if (!length(gene_reactions)) {
        stop(
          paste(
            "No complete GPR alternative is covered by `core_genes`;",
            "use `gene_match = \"any_direct\"` only when partial-complex",
            "selection is intentional."
          ),
          call. = FALSE
        )
      }
    }
    gene_source <- lapply(gene_reactions, function(reaction) {
      sort(unique(gpr$gene[
        gpr$reaction_id == reaction & gpr$gene %in% requested_genes
      ]))
    })
    names(gene_source) <- gene_reactions
  }

  reactions <- union(requested_reactions, gene_reactions)
  source <- vapply(reactions, function(reaction) {
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
    if (is.null(genes) || !length(genes)) NA_character_ else paste(genes, collapse = ";")
  }, character(1))

  data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    gene = mapped_genes,
    reaction_id = reactions,
    is_core = TRUE,
    selection_source = source,
    stringsAsFactors = FALSE
  )
}

.rc_build_target_union_definition <- function(
    gem, core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L) {
  gene_match <- match.arg(gene_match)
  expansion_mode <- match.arg(expansion_mode)
  core <- .rc_target_union_core_rows(
    gem = gem,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
  expanded <- rc_expand_meta_module_reactions(
    gem = gem,
    core_reactions = core,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )
  membership <- expanded$reaction_membership
  membership$score_target <- membership$is_core %in% TRUE
  membership$model_only <- !membership$score_target
  core <- core[
    match(unique(as.character(membership$reaction_id[membership$score_target])),
          core$reaction_id),
    , drop = FALSE
  ]
  rownames(core) <- NULL
  rownames(membership) <- NULL
  summary <- expanded$summary
  summary$n_score_targets <- nrow(core)
  summary$n_model_only_reactions <- sum(membership$model_only)
  summary$gene_match <- gene_match
  summary$expansion_mode <- expansion_mode

  list(
    global_core_reactions = core,
    global_reaction_membership = membership,
    summary = summary,
    crossref_maps = expanded$crossref_maps,
    params = list(
      gene_match = gene_match,
      expansion_mode = expansion_mode,
      max_iterations = as.integer(max_iterations),
      score_targets = unique(as.character(core$reaction_id)),
      model_only_reactions = unique(as.character(
        membership$reaction_id[membership$model_only]
      ))
    )
  )
}

#' Score manually selected core reactions in an annotation-expanded union GEM
#'
#' This optional step accepts one or more reaction IDs and/or genes that directly
#' map to GEM GPR rules. Selected reactions are the only directional scoring
#' targets. Reactions sharing a core subsystem, KEGG or Reactome reaction ID, or
#' master Rhea ID are added to one union GEM as model-only reactions. Existing
#' add-only FASTCORE completion may add further support reactions when required
#' for parent-GEM feasibility.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()] or another compatible
#'   Layer 1 object containing reaction expression and unit metadata.
#' @param gem Validated GEM with reaction metadata and a parsed GPR table.
#' @param medium_scenarios Shared medium table used for all scored units.
#' @param outdir Persistent output directory.
#' @param core_reaction_ids Character vector of reactions to score.
#' @param core_genes Character vector of genes used to select directly associated
#'   reactions.
#' @param gene_match `"complete_gpr"` requires all genes in at least one GPR AND
#'   group. `"any_direct"` permits intentionally selecting a reaction from one
#'   directly mapped subunit.
#' @param subsystem_table Optional external reaction-to-subsystem table.
#' @param expansion_mode Annotation expansion mode passed to
#'   [rc_expand_meta_module_reactions()].
#' @param max_iterations Maximum fixed-point annotation-expansion iterations.
#' @param condition_col,sample_col,celltype_col Layer 1 metadata columns.
#' @param layer2_args Additional arguments passed to [rc_run_microcompass()].
#' @param parallel Whether LP tasks may run in parallel.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @return A `regcompass_target_union_step` containing the selected cores,
#'   model-only union membership, expansion summary, and microCOMPASS result.
#' @export
rc_regcompass_step_target_union <- function(
    layer1, gem, medium_scenarios, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L,
    condition_col = "condition", sample_col = "sample_id",
    celltype_col = "cell_type",
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {
  gene_match <- match.arg(gene_match)
  expansion_mode <- match.arg(expansion_mode)
  if (!is.list(layer1) || is.null(layer1$reaction_expression) ||
      is.null(layer1$unit_meta)) {
    stop(
      "`layer1` must contain reaction_expression and unit_meta.",
      call. = FALSE
    )
  }
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  reserved <- intersect(names(layer2_args), c(
    "layer1", "gem", "target_reactions", "medium_scenarios", "mode",
    "reaction_membership", "core_reactions", "unit", "condition_col",
    "sample_col", "celltype_col", "parallel", "BPPARAM"
  ))
  if (length(reserved)) {
    stop(
      "`layer2_args` cannot override target-union workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }

  definition <- .rc_build_target_union_definition(
    gem = gem,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  layer2_args$model_params$cache_dir <- file.path(outdir, "model_cache")
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = definition$global_core_reactions$reaction_id,
    medium_scenarios = medium_scenarios,
    mode = "meta_module_gem",
    reaction_membership = definition$global_reaction_membership,
    core_reactions = definition$global_core_reactions,
    unit = "metacell",
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  scored <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(w) {
      if (grepl(
        "Metacell-level scores are descriptive pseudo-observations",
        conditionMessage(w), fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  scored$params$target_scope <- "manual_core_annotation_expanded_union_gem"
  scored$params$n_score_targets <- nrow(definition$global_core_reactions)
  scored$params$n_model_only_reactions <- sum(
    definition$global_reaction_membership$model_only
  )

  .rc_mm_write_tsv_gz(
    definition$global_core_reactions,
    file.path(outdir, "selected_core_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$global_reaction_membership,
    file.path(outdir, "union_reaction_membership.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$summary,
    file.path(outdir, "union_summary.tsv.gz")
  )
  rc_export_microcompass(scored, file.path(outdir, "scores"))
  answer <- c(definition, list(microcompass = scored))
  class(answer) <- c("regcompass_target_union_step", "list")
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
