#' Map GRN metabolic gene modules to core Human-GEM reactions
#' @export
rc_map_meta_module_core_reactions <- function(gene_nodes, gpr_table) {
  if (!is.data.frame(gene_nodes) || !all(c("sample_id", "gene", "module_id") %in% colnames(gene_nodes))) {
    stop("`gene_nodes` must contain sample_id, gene and module_id.", call. = FALSE)
  }
  if (!is.data.frame(gpr_table) || !all(c("reaction_id", "gene") %in% colnames(gpr_table))) {
    stop("`gpr_table` must contain reaction_id and gene.", call. = FALSE)
  }
  gpr <- gpr_table
  gpr$gene <- toupper(as.character(gpr$gene))
  nodes <- gene_nodes
  nodes$gene <- toupper(as.character(nodes$gene))
  out <- merge(nodes, unique(gpr[, c("reaction_id", "gene"), drop = FALSE]), by = "gene", all = FALSE, sort = FALSE)
  out$is_core <- TRUE
  out$inclusion_stage <- "core_grn_gene"
  out <- unique(out[, c("sample_id", "module_id", "gene", "reaction_id", "is_core", "inclusion_stage"), drop = FALSE])
  rownames(out) <- NULL
  out
}

#' Expand core reactions into GRN-defined reaction meta-modules
#'
#' Expansion is deliberately ordered and non-recursive by default: core reaction
#' subsystems, then shared KEGG/Reactome identifiers, then shared master-Rhea IDs.
#' @export
rc_expand_meta_module_reactions <- function(gem,
                                            core_reactions,
                                            subsystem_table = NULL,
                                            expansion_mode = c("ordered_once", "fixed_point"),
                                            max_iterations = 10L) {
  expansion_mode <- match.arg(expansion_mode)
  required <- c("sample_id", "module_id", "gene", "reaction_id")
  if (!is.data.frame(core_reactions) || !all(required %in% colnames(core_reactions))) {
    stop("`core_reactions` must contain sample_id, module_id, gene and reaction_id.", call. = FALSE)
  }
  maps <- rc_reaction_crossref_maps(gem, subsystem_table = subsystem_table)
  if (!nrow(maps$subsystem)) {
    stop("No usable reaction-to-subsystem annotations were found. Use `rc_prepare_human2_gem_v12()` or supply `subsystem_table`.", call. = FALSE)
  }
  valid_reactions <- colnames(rc_validate_gem(gem)$S)
  core_reactions <- core_reactions[core_reactions$reaction_id %in% valid_reactions, , drop = FALSE]
  if (!nrow(core_reactions)) stop("No GRN genes mapped to valid GEM reactions.", call. = FALSE)
  groups <- split(seq_len(nrow(core_reactions)), paste(core_reactions$sample_id, core_reactions$module_id, sep = "\001"))
  membership_rows <- list()
  summary_rows <- list()

  for (key in names(groups)) {
    core <- unique(as.character(core_reactions$reaction_id[groups[[key]]]))
    sample <- as.character(core_reactions$sample_id[groups[[key]][[1L]]])
    module <- as.character(core_reactions$module_id[groups[[key]][[1L]]])
    included <- core
    reasons <- stats::setNames(rep("core_grn_gene", length(core)), core)
    source_ids <- stats::setNames(rep(NA_character_, length(core)), core)

    add_reactions <- function(rxns, reason, source_map = NULL) {
      rxns <- intersect(.rc_mm_trim_unique(rxns), valid_reactions)
      new <- setdiff(rxns, included)
      if (length(new)) {
        included <<- c(included, new)
        reasons[new] <<- reason
        if (!is.null(source_map)) source_ids[new] <<- source_map[new]
      }
      invisible(new)
    }

    expand_once <- function() {
      before <- length(included)
      reactions_in_subsystems <- function(subsystems) {
        unique(maps$subsystem$reaction_id[maps$subsystem$subsystem_id %in% subsystems])
      }
      subsystems_for_reactions <- function(rxns) {
        unique(maps$subsystem$subsystem_id[maps$subsystem$reaction_id %in% rxns])
      }
      source_for_subsystem_members <- function(subsystems, labels_by_subsystem) {
        rxns <- reactions_in_subsystems(subsystems)
        stats::setNames(vapply(rxns, function(r) {
          ss <- unique(maps$subsystem$subsystem_id[
            maps$subsystem$reaction_id == r & maps$subsystem$subsystem_id %in% subsystems])
          labels <- unique(unlist(labels_by_subsystem[ss], use.names = FALSE))
          paste(labels[!is.na(labels) & nzchar(labels)], collapse = ";")
        }, character(1)), rxns)
      }

      core_subsystems <- subsystems_for_reactions(core)
      subsystem_rxns <- reactions_in_subsystems(core_subsystems)
      subsystem_labels <- stats::setNames(lapply(core_subsystems, function(ss) paste0("SUBSYSTEM:", ss)), core_subsystems)
      subsystem_source <- source_for_subsystem_members(core_subsystems, subsystem_labels)
      add_reactions(subsystem_rxns, "same_core_subsystem", subsystem_source)

      current_for_db <- included
      kegg_ids <- unique(maps$kegg$kegg_id[maps$kegg$reaction_id %in% current_for_db])
      reactome_ids <- unique(maps$reactome$reactome_id[maps$reactome$reaction_id %in% current_for_db])
      all_subsystems <- unique(maps$subsystem$subsystem_id)
      db_labels <- lapply(all_subsystems, function(ss) {
        ss_rxns <- reactions_in_subsystems(ss)
        kid <- intersect(unique(maps$kegg$kegg_id[maps$kegg$reaction_id %in% ss_rxns]), kegg_ids)
        rid <- intersect(unique(maps$reactome$reactome_id[maps$reactome$reaction_id %in% ss_rxns]), reactome_ids)
        unique(c(paste0("SUBSYSTEM:", ss), paste0("KEGG:", kid), paste0("REACTOME:", rid)))
      })
      names(db_labels) <- all_subsystems
      db_subsystems <- all_subsystems[vapply(db_labels, function(z) any(grepl("^(KEGG|REACTOME):", z)), logical(1))]
      db_rxns <- reactions_in_subsystems(db_subsystems)
      db_source <- source_for_subsystem_members(db_subsystems, db_labels)
      add_reactions(db_rxns, "shared_kegg_or_reactome_subsystem", db_source)

      current_for_rhea <- included
      master_ids <- unique(maps$rhea_master$rhea_master_id[maps$rhea_master$reaction_id %in% current_for_rhea])
      rhea_anchor_rxns <- unique(maps$rhea_master$reaction_id[maps$rhea_master$rhea_master_id %in% master_ids])
      rhea_subsystems <- subsystems_for_reactions(rhea_anchor_rxns)
      rhea_labels <- lapply(rhea_subsystems, function(ss) {
        ss_rxns <- reactions_in_subsystems(ss)
        mids <- intersect(unique(maps$rhea_master$rhea_master_id[maps$rhea_master$reaction_id %in% ss_rxns]), master_ids)
        unique(c(paste0("SUBSYSTEM:", ss), paste0("RHEA_MASTER:", mids)))
      })
      names(rhea_labels) <- rhea_subsystems
      rhea_rxns <- reactions_in_subsystems(rhea_subsystems)
      rhea_source <- source_for_subsystem_members(rhea_subsystems, rhea_labels)
      add_reactions(rhea_rxns, "shared_master_rhea_subsystem", rhea_source)
      length(included) > before
    }

    changed <- expand_once()
    iteration <- 1L
    if (identical(expansion_mode, "fixed_point")) {
      while (isTRUE(changed) && iteration < as.integer(max_iterations)) {
        iteration <- iteration + 1L
        changed <- expand_once()
      }
    }
    member <- data.frame(sample_id = sample, module_id = module,
                         reaction_id = included, is_core = included %in% core,
                         inclusion_stage = unname(reasons[included]),
                         source_annotation = unname(source_ids[included]),
                         expansion_mode = expansion_mode,
                         stringsAsFactors = FALSE)
    membership_rows[[key]] <- member
    summary_rows[[key]] <- data.frame(
      sample_id = sample, module_id = module,
      n_core_genes = length(unique(core_reactions$gene[groups[[key]]])),
      n_core_reactions = length(core), n_reactions = length(included),
      n_subsystem_added = sum(member$inclusion_stage == "same_core_subsystem"),
      n_database_added = sum(member$inclusion_stage == "shared_kegg_or_reactome_subsystem"),
      n_rhea_added = sum(member$inclusion_stage == "shared_master_rhea_subsystem"),
      iterations = iteration, stringsAsFactors = FALSE)
  }
  list(reaction_membership = do.call(rbind, membership_rows),
       summary = do.call(rbind, summary_rows),
       crossref_maps = maps)
}

#' Build a module-meso-GEM from one GRN-defined reaction meta-module
#' @export
rc_build_meta_module_gem <- function(gem,
                                     reaction_membership,
                                     sample_id,
                                     module_id,
                                     medium_table = NULL,
                                     condition = NULL,
                                     include_one_hop = FALSE,
                                     include_transport = TRUE,
                                     include_exchange = TRUE,
                                     include_protected = TRUE,
                                     currency_metabolites = NULL,
                                     max_reactions = 3000,
                                     strict_closure = FALSE) {
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(reaction_membership) || !all(required %in% colnames(reaction_membership))) {
    stop("`reaction_membership` must contain sample_id, module_id and reaction_id.", call. = FALSE)
  }
  biological <- unique(as.character(reaction_membership$reaction_id[
    as.character(reaction_membership$sample_id) == as.character(sample_id) &
      as.character(reaction_membership$module_id) == as.character(module_id)]))
  if (!length(biological)) stop("No reactions found for the requested sample/module.", call. = FALSE)
  gv <- rc_validate_gem(gem)
  gem2 <- gem
  meta <- gem2$reaction_meta
  if (is.null(meta)) meta <- data.frame(reaction_id = gv$reactions, stringsAsFactors = FALSE)
  meta <- meta[match(gv$reactions, as.character(meta$reaction_id)), , drop = FALSE]
  module_key <- as.character(module_id)
  meta$grn_meta_module <- ifelse(meta$reaction_id %in% biological, module_key, "OUTSIDE_GRN_META_MODULE")
  gem2$reaction_meta <- meta
  out <- rc_build_module_meso_gem(
    gem = gem2, module_id = module_key, medium_table = medium_table,
    condition = condition, module_col = "grn_meta_module",
    include_one_hop = include_one_hop, include_transport = include_transport,
    include_exchange = include_exchange, include_protected = include_protected,
    currency_metabolites = currency_metabolites, max_reactions = max_reactions,
    strict_closure = strict_closure
  )
  out$sample_id <- as.character(sample_id)
  out$grn_module_id <- as.character(module_id)
  out$reaction_meta$biological_meta_module_member <- out$reaction_meta$reaction_id %in% biological
  out$reaction_meta$support_only <- !out$reaction_meta$biological_meta_module_member
  out
}

#' Load the retained metacell Seurat object from a formal RegCompass run
#' @export
rc_load_metacell_object_from_run <- function(run_dir,
                                             retained_metacell_ids = NULL,
                                             rna_assay = "RNA",
                                             atac_assay = "ATAC") {
  metacell_root <- file.path(run_dir, "01_metacells")
  object_files <- list.files(metacell_root, pattern = "metacell_object\\.rds$", recursive = TRUE, full.names = TRUE)
  if (!length(object_files)) stop("No metacell_object.rds files were found under `01_metacells`.", call. = FALSE)
  objects <- lapply(object_files, readRDS)
  object_cells <- unlist(lapply(objects, colnames), use.names = FALSE)
  if (anyDuplicated(object_cells)) stop("Saved metacell objects contain duplicated metacell IDs.", call. = FALSE)
  if (!is.null(retained_metacell_ids)) {
    retained_metacell_ids <- .rc_mm_trim_unique(retained_metacell_ids)
    missing <- setdiff(retained_metacell_ids, object_cells)
    if (length(missing)) stop("Retained metacell IDs are missing from saved objects: ", paste(utils::head(missing, 10L), collapse = ", "), call. = FALSE)
    objects <- lapply(objects, function(x) {
      keep <- intersect(colnames(x), retained_metacell_ids)
      if (!length(keep)) return(NULL)
      subset(x, cells = keep)
    })
    objects <- objects[!vapply(objects, is.null, logical(1))]
  }
  objects <- lapply(objects, .rc_clear_signac_fragments, atac_assay = atac_assay)
  out <- if (length(objects) == 1L) objects[[1L]] else Reduce(function(a, b) merge(x = a, y = b, merge.data = FALSE), objects)
  if (!is.null(retained_metacell_ids)) out <- subset(out, cells = retained_metacell_ids)
  out
}


#' Run the integrated RegCompassR v1.2 workflow
#' @export
rc_run_regcompass_v12 <- function(object,
                                  gem,
                                  outdir,
                                  pfm,
                                  genome,
                                  fragment_files = NULL,
                                  sample_col = "sample_id",
                                  condition_col = "condition",
                                  celltype_col = "cell_type",
                                  rna_assay = "RNA",
                                  atac_assay = "ATAC",
                                  metacell_args = list(),
                                  pando_args = list()) {
  if (is.null(gem$gpr_table)) stop("`gem` must contain `gpr_table`.", call. = FALSE)
  formal_defaults <- list(object = object, gpr_table = gem$gpr_table, outdir = outdir,
                          fragment_files = fragment_files, sample_col = sample_col,
                          condition_col = condition_col, celltype_col = celltype_col,
                          rna_assay = rna_assay, atac_assay = atac_assay)
  formal_defaults[names(metacell_args)] <- NULL
  layer1 <- do.call(rc_run_regcompass_multiome_metacell, c(formal_defaults, metacell_args))
  retained <- as.character(layer1$metacell_meta$metacell_id)
  mc_object <- rc_load_metacell_object_from_run(outdir, retained_metacell_ids = retained,
                                                 rna_assay = rna_assay, atac_assay = atac_assay)
  single_cell_genes <- rownames(.rc_get_assay_counts(object, rna_assay))
  pando_defaults <- list(metacell_object = mc_object, gem = gem,
                         outdir = file.path(outdir, "04_pando_meta_modules"),
                         pfm = pfm, genome = genome, sample_col = sample_col,
                         single_cell_genes = single_cell_genes,
                         rna_assay = rna_assay, atac_assay = atac_assay)
  pando_defaults[names(pando_args)] <- NULL
  layer1$grn_meta_modules <- do.call(rc_run_pando_meta_modules, c(pando_defaults, pando_args))
  layer1$schema_version <- "regcompass_v1.2"
  saveRDS(layer1, file.path(outdir, "regcompass_v1.2_result.rds"))
  layer1
}
