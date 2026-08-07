# Merge biological meta-modules within each cell type only.

.rc_merge_meta_modules_by_cell_type <- function(
    condition_modules,
    celltype_col = "cell_type",
    condition_col = "condition") {
  if (!is.list(condition_modules)) {
    stop("`condition_modules` must be a list.", call. = FALSE)
  }
  valid_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (!valid_name(celltype_col) || !valid_name(condition_col)) {
    stop("Cell-type and condition column names must be non-empty strings.",
         call. = FALSE)
  }

  names_to_merge <- c(
    "supported_metabolic_genes", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  tables <- lapply(names_to_merge, function(name) {
    value <- condition_modules[[name]]
    if (is.data.frame(value)) value else data.frame()
  })
  names(tables) <- names_to_merge

  normalize_scope <- function(tab, label) {
    if (!is.data.frame(tab) || !nrow(tab) ||
        !all(c(celltype_col, "reaction_id") %in% colnames(tab))) {
      stop(
        "`", label, "` must contain non-empty `", celltype_col,
        "` and `reaction_id` columns.",
        call. = FALSE
      )
    }
    tab[[celltype_col]] <- trimws(as.character(tab[[celltype_col]]))
    tab$reaction_id <- trimws(as.character(tab$reaction_id))
    if (anyNA(tab[[celltype_col]]) || any(!nzchar(tab[[celltype_col]])) ||
        anyNA(tab$reaction_id) || any(!nzchar(tab$reaction_id))) {
      stop("Meta-module cell types and reaction IDs must be complete.",
           call. = FALSE)
    }
    tab
  }

  core <- normalize_scope(tables$core_gene_reaction, "core_gene_reaction")
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  biological <- normalize_scope(
    tables$reaction_membership, "reaction_membership"
  )
  cell_types <- sort(unique(as.character(biological[[celltype_col]])))
  core_cell_types <- unique(as.character(core[[celltype_col]]))
  missing_core <- setdiff(cell_types, core_cell_types)
  if (length(missing_core)) {
    stop(
      "No complete-GPR core reactions were produced for cell types: ",
      paste(missing_core, collapse = ", "),
      call. = FALSE
    )
  }

  subset_table <- function(tab, cell_type) {
    if (!is.data.frame(tab) || !nrow(tab)) return(tab)
    if (!celltype_col %in% colnames(tab)) return(tab[0, , drop = FALSE])
    value <- trimws(as.character(tab[[celltype_col]]))
    tab[value == cell_type, , drop = FALSE]
  }

  build_one <- function(cell_type) {
    one <- lapply(tables, subset_table, cell_type = cell_type)
    one_core <- core[core[[celltype_col]] == cell_type, , drop = FALSE]
    one_biological <- biological[
      biological[[celltype_col]] == cell_type, , drop = FALSE
    ]
    core_ids <- unique(as.character(one_core$reaction_id))
    biological_ids <- unique(as.character(one_biological$reaction_id))
    if (!length(core_ids) || !length(biological_ids)) {
      stop("Cell type `", cell_type, "` has an empty meta-module catalogue.",
           call. = FALSE)
    }
    missing <- setdiff(core_ids, biological_ids)
    if (length(missing)) {
      stop(
        "Cell type `", cell_type,
        "` has core reactions missing from its biological membership: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    catalogue_id <- paste0(
      "CELLTYPE_META_MODULES::",
      utils::URLencode(cell_type, reserved = TRUE)
    )
    scoped_core <- data.frame(
      catalogue_id = catalogue_id,
      cell_type = cell_type,
      reaction_id = core_ids,
      is_core = TRUE,
      stringsAsFactors = FALSE
    )
    names(scoped_core)[names(scoped_core) == "cell_type"] <- celltype_col
    scoped_membership <- data.frame(
      catalogue_id = catalogue_id,
      cell_type = cell_type,
      reaction_id = biological_ids,
      is_core = biological_ids %in% core_ids,
      inclusion_stage = ifelse(
        biological_ids %in% core_ids,
        "celltype_union_core",
        "celltype_union_biological_member"
      ),
      stringsAsFactors = FALSE
    )
    names(scoped_membership)[names(scoped_membership) == "cell_type"] <-
      celltype_col

    source_groups <- if ("group_id" %in% colnames(one_biological)) {
      value <- trimws(as.character(one_biological$group_id))
      sort(unique(value[!is.na(value) & nzchar(value)]))
    } else {
      character()
    }
    source_conditions <- if (condition_col %in% colnames(one_biological)) {
      value <- trimws(as.character(one_biological[[condition_col]]))
      sort(unique(value[!is.na(value) & nzchar(value)]))
    } else {
      character()
    }

    one$biological_reaction_membership <- one_biological
    one$merged_core_reactions <- scoped_core
    one$merged_reaction_membership <- scoped_membership
    one$schema_version <- "regcompass_celltype_meta_modules_v1"
    one$cell_type <- cell_type
    one$celltype_col <- celltype_col
    one$condition_col <- condition_col
    one$source_group_ids <- source_groups
    one$source_conditions <- source_conditions
    one$core_definition <-
      "condition_specific_complete_gpr_cores_unioned_within_cell_type"
    one$merge_source <-
      "condition_meta_module_reactions_deduplicated_within_cell_type"
    one$merge_scope <- "cell_type"
    one$is_gem <- FALSE
    one$fastcore_applied <- FALSE
    one
  }

  catalogues <- stats::setNames(lapply(cell_types, build_one), cell_types)
  bind <- function(field) {
    .rc_bind_frames_fill(lapply(catalogues, `[[`, field))
  }
  out <- tables
  out$biological_reaction_membership <- biological
  out$cell_type_catalogues <- catalogues
  out$merged_core_reactions <- bind("merged_core_reactions")
  out$merged_reaction_membership <- bind("merged_reaction_membership")
  out$schema_version <- "regcompass_celltype_merged_meta_modules_v1"
  out$celltype_col <- celltype_col
  out$condition_col <- condition_col
  out$cell_types <- cell_types
  out$core_definition <-
    "condition_specific_complete_gpr_cores_unioned_within_cell_type"
  out$merge_source <-
    "condition_meta_module_reactions_deduplicated_within_cell_type_only"
  out$merge_scope <- "cell_type"
  out$cross_celltype_merge <- FALSE
  out$is_gem <- FALSE
  out$fastcore_applied <- FALSE
  out
}
