from pathlib import Path


def insert_once(path, marker, insertion):
    text = path.read_text()
    if text.count(marker) != 1:
        raise RuntimeError(f'{path}: marker count is {text.count(marker)}: {marker}')
    path.write_text(text.replace(marker, insertion + marker, 1))


def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise RuntimeError(f'{path}: expected one match, found {text.count(old)}')
    path.write_text(text.replace(old, new, 1))


target = Path('R/target_union.R')
helper = r'''.rc_validate_target_union_anchor_request <- function(
    gem, merged_core_reactions,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct")) {
  gene_match <- match.arg(gene_match)
  if (!is.data.frame(merged_core_reactions) ||
      !"reaction_id" %in% colnames(merged_core_reactions)) {
    stop("`merged_core_reactions` must contain reaction_id.",
         call. = FALSE)
  }
  available_core <- .rc_target_union_normalize_ids(
    merged_core_reactions$reaction_id
  )
  if (!length(available_core)) {
    stop("No original cell-type core reaction is available.", call. = FALSE)
  }
  requested_reactions <- .rc_target_union_normalize_ids(core_reaction_ids)
  requested_genes <- toupper(.rc_target_union_normalize_ids(core_genes))
  if (!length(requested_reactions) && !length(requested_genes)) {
    stop(
      "Supply at least one `core_reaction_ids` or `core_genes` value.",
      call. = FALSE
    )
  }

  if (length(requested_reactions)) {
    .rc_target_union_core_rows(
      gem = gem,
      available_core_reactions = available_core,
      core_reaction_ids = requested_reactions,
      gene_match = gene_match
    )
  }
  if (length(requested_genes)) {
    .rc_target_union_core_rows(
      gem = gem,
      available_core_reactions = available_core,
      core_genes = requested_genes,
      gene_match = gene_match
    )
  }
  invisible(list(
    core_reaction_ids = requested_reactions,
    core_genes = requested_genes,
    gene_match = gene_match
  ))
}

'''
insert_once(
    target,
    '.rc_build_target_union_definition <- function(',
    helper
)

validation_marker = '''  if (!is.list(cached_reaction_ids) || is.null(names(cached_reaction_ids))) {
    stop("Cached reaction IDs must be a named list by cell type.",
         call. = FALSE)
  }
'''
validation_replacement = validation_marker + '''  .rc_validate_target_union_anchor_request(
    gem = gem,
    merged_core_reactions = merged_core_reactions,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
'''
replace_once(target, validation_marker, validation_replacement)

old_roxygen = '''#' Score directly database-linked non-core reactions in final union GEMs
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
'''
new_roxygen = '''#' Score directly linked non-core reactions in cell-type union GEMs
#'
#' Selected reactions are mapping anchors only. Direct KEGG-, Reactome-, or
#' master-Rhea-linked non-core reactions are scored by reusing the exact Stage 5
#' union GEMs for the corresponding cell types and media. Candidate availability
#' is intersected across media within each cell type. Models from different cell
#' types are never merged, and FASTCORE is not rerun.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer2 Output from [rc_regcompass_step_layer2()] with
#'   `model_mode = "meta_module_gem"`.
#' @param gem The same GEM used for the original run.
#' @param outdir Output directory.
#' @param core_reaction_ids GEM reaction IDs used as direct mapping anchors.
#'   The historical argument name is retained for compatibility. Every supplied
#'   ID must exist in `gem`; it is evaluated only in cell types where it is
#'   structurally available.
#' @param core_genes Genes used to resolve original core anchors through GPRs
#'   within each cell type. Every supplied gene selector must resolve to at least
#'   one original cell-type core under `gene_match`.
#' @param gene_match Require a complete GPR group or allow any direct gene match.
#' @param layer2_args Optional `omega`, `target_direction`, `solver`, and
#'   `flux_threshold` overrides. Scoring LPs have no construction time limit.
#' @param parallel Whether to parallelize reused cell-type-model by matching-
#'   metacell tasks.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param progress Whether to display stage progress.
#' @return A `regcompass_target_union_step` containing cell-type-scoped anchors,
#'   direct database-linked non-core targets, exact model provenance, and LP
#'   scores.
#' @export
'''
replace_once(target, old_roxygen, new_roxygen)

old_counts = '''  scored$params$target_scope <-
    "direct_kegg_reactome_master_rhea_noncore_only"
  scored$params$n_selected_core <- nrow(definition$selected_core_reactions)
  scored$params$n_merged_core_reactions_not_rescored <-
    length(definition$params$merged_core_reactions_not_rescored)
  scored$params$n_cached_union_unavailable_reactions <-
    length(definition$params$cached_union_unavailable_reactions)
  scored$params$n_expanded_score_targets <-
    nrow(definition$params$score_targets)
'''
new_counts = '''  scored$params$target_scope <-
    "direct_database_crossrefs_within_cell_type_only"
  scored$params$n_selected_core <- nrow(definition$selected_core_reactions)
  relation_catalogue <- definition$expanded_reaction_catalog
  scored$params$n_merged_core_reactions_not_rescored <- length(unique(
    as.character(relation_catalogue$reaction_id[
      relation_catalogue$merged_catalogue_is_core %in% TRUE
    ])
  ))
  scored$params$n_cached_union_unavailable_reactions <- length(unique(
    as.character(relation_catalogue$reaction_id[
      !relation_catalogue$available_in_all_cached_union_gems %in% TRUE
    ])
  ))
  scored$params$n_expanded_score_targets <-
    nrow(definition$params$score_targets)
'''
replace_once(target, old_counts, new_counts)

step2 = Path('R/step_layer2.R')
old_layer2_roxygen = '''#' Build medium-specific structural models and run directional LP scoring
#'
#' The primary route uses the current fixed-dictionary regulatory evidence.
#' Historical `_oof`, `common`, and `condition_unique` fields are retained as
#' compatibility aliases. RNA-only scoring remains an interpretation control,
#' and all routes reuse the exact same medium-specific structural model.
#'
#' @export
'''
new_layer2_roxygen = '''#' Build cell-type and medium structural models and score directional LPs
#'
#' With `model_mode = "meta_module_gem"`, conditions are unioned only within the
#' same cell type. One union GEM and one independent FASTCORE completion are
#' created for every cell-type and medium combination. Conditions and metacells
#' share a model only when their cell type matches. Historical `_oof`, `common`,
#' and `condition_unique` fields remain compatibility aliases; RNA-only scoring
#' is an interpretation control. The `full_gem` route uses a separate engine.
#'
#' @export
'''
replace_once(step2, old_layer2_roxygen, new_layer2_roxygen)

test = Path('tests/testthat/test-target-union.R')
test_text = test.read_text()
marker = '''  expect_error(
    .rc_target_union_core_rows(
      gem,
      available_core_reactions = available,
      core_reaction_ids = "missing"
    ),
    "absent from the GEM"
  )
'''
addition = marker + '''  expect_error(
    .rc_build_target_union_definition(
      gem = gem,
      merged_core_reactions = target_union_merged_core(),
      merged_reaction_membership = target_union_merged_membership(),
      core_reaction_ids = "missing",
      cached_reaction_ids = list(C = paste0("R", 1:7)),
      celltype_col = "cell_type"
    ),
    "absent from the GEM"
  )
  expect_error(
    .rc_build_target_union_definition(
      gem = gem,
      merged_core_reactions = target_union_merged_core(),
      merged_reaction_membership = target_union_merged_membership(),
      core_genes = "NOT_A_GPR_GENE",
      cached_reaction_ids = list(C = paste0("R", 1:7)),
      celltype_col = "cell_type"
    ),
    "do not map to GEM GPR rules"
  )
'''
if test_text.count(marker) != 1:
    raise RuntimeError('target validation test marker was not found once')
test.write_text(test_text.replace(marker, addition, 1))
