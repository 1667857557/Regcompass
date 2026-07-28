.rc_summarize_supported_metabolic_genes <- function(
    grn_result, metabolic_genes) {
  if (!is.list(grn_result) ||
      !is.data.frame(grn_result$tf_peak_gene_significant)) {
    stop("`grn_result` is not a valid single-cell GRN result.", call. = FALSE)
  }
  group_cols <- as.character(grn_result$group_cols)
  required <- c("group_id", group_cols, "tf", "target", "region")
  active <- grn_result$tf_peak_gene_significant
  missing <- setdiff(required, colnames(active))
  if (length(missing)) {
    stop(
      "Active Pando edge table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  metabolic_genes <- unique(toupper(.rc_mm_trim_unique(metabolic_genes)))
  active$target <- toupper(trimws(as.character(active$target)))
  active$tf <- toupper(trimws(as.character(active$tf)))
  active$region <- trimws(as.character(active$region))
  active <- active[
    active$target %in% metabolic_genes,
    , drop = FALSE
  ]
  if (!nrow(active)) {
    stop(
      paste(
        "No Human-GEM metabolic target genes have active Pando",
        "TF-peak-gene evidence."
      ),
      call. = FALSE
    )
  }

  key <- paste(active$group_id, active$target, sep = "\001")
  rows <- split(seq_len(nrow(active)), key)
  summary_rows <- lapply(rows, function(index) {
    one <- active[index, , drop = FALSE]
    group_id <- as.character(one$group_id[[1L]])
    target <- as.character(one$target[[1L]])
    estimate <- if ("estimate" %in% colnames(one)) {
      suppressWarnings(as.numeric(one$estimate))
    } else {
      rep(NA_real_, nrow(one))
    }
    padj <- if ("padj" %in% colnames(one)) {
      suppressWarnings(as.numeric(one$padj))
    } else {
      rep(NA_real_, nrow(one))
    }
    rsq <- if ("rsq" %in% colnames(one)) {
      suppressWarnings(as.numeric(one$rsq))
    } else {
      rep(NA_real_, nrow(one))
    }
    finite_min <- function(value) {
      value <- value[is.finite(value)]
      if (length(value)) min(value) else NA_real_
    }
    finite_max <- function(value) {
      value <- value[is.finite(value)]
      if (length(value)) max(value) else NA_real_
    }
    group_values <- one[1L, group_cols, drop = FALSE]
    data.frame(
      group_id = group_id,
      group_values,
      sample_id = group_id,
      module_id = paste0(group_id, "::SUPPORTED_METABOLIC_GENES"),
      gene = target,
      n_active_edges = nrow(one),
      # Retained as a compatibility alias for RegCompassR <= 1.9.0.
      n_significant_edges = nrow(one),
      n_regulating_tfs = length(unique(one$tf[nzchar(one$tf)])),
      n_regulatory_regions = length(unique(one$region[nzchar(one$region)])),
      min_padj = finite_min(padj),
      max_abs_estimate = finite_max(abs(estimate)),
      max_model_rsq = finite_max(rsq),
      n_positive_edges = sum(is.finite(estimate) & estimate > 0),
      n_negative_edges = sum(is.finite(estimate) & estimate < 0),
      evidence_definition = "active_pando_tf_peak_gene_target",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  answer <- unique(do.call(rbind, summary_rows))
  rownames(answer) <- NULL
  answer
}

.rc_build_condition_meta_modules <- function(
    grn_result, gem, outdir, meta_module_args = list()) {
  if (!is.list(meta_module_args)) {
    stop("`meta_module_args` must be a list.", call. = FALSE)
  }
  allowed <- "subsystem_table"
  unknown <- setdiff(names(meta_module_args), allowed)
  if (length(unknown)) {
    stop(
      "Unknown `meta_module_args` fields: ",
      paste(unknown, collapse = ", "),
      ". Allowed field: `subsystem_table`.",
      call. = FALSE
    )
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_cols <- grn_result$group_cols
  display_cols <- c("group_id", group_cols)
  module_cols <- unique(c(display_cols, "sample_id", "module_id"))
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)

  supported <- .rc_summarize_supported_metabolic_genes(
    grn_result = grn_result,
    metabolic_genes = metabolic_genes
  )
  core <- rc_map_meta_module_core_reactions(
    supported[, c("sample_id", "module_id", "gene"), drop = FALSE],
    gem$gpr_table
  )
  if (nrow(core)) {
    core <- merge(
      core,
      supported,
      by = c("sample_id", "module_id", "gene"),
      all.x = TRUE,
      sort = FALSE
    )
    core <- core[, c(
      display_cols,
      setdiff(colnames(core), display_cols)
    ), drop = FALSE]
  }
  if (!nrow(core) || !any(core$is_core %in% TRUE)) {
    stop(
      paste(
        "Significant Pando-supported metabolic genes did not completely",
        "satisfy any Human-GEM GPR branch."
      ),
      call. = FALSE
    )
  }

  expanded <- rc_expand_meta_module_reactions(
    gem,
    core,
    subsystem_table = meta_module_args$subsystem_table %||% NULL
  )
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(
      expanded$reaction_membership,
      unique(supported[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    expanded$reaction_membership <- expanded$reaction_membership[, c(
      display_cols,
      setdiff(colnames(expanded$reaction_membership), display_cols)
    ), drop = FALSE]
  }
  if (nrow(expanded$summary)) {
    expanded$summary <- merge(
      expanded$summary,
      unique(supported[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    expanded$summary <- expanded$summary[, c(
      display_cols,
      setdiff(colnames(expanded$summary), display_cols)
    ), drop = FALSE]
  }

  out <- c(grn_result, list(
    supported_metabolic_genes = supported,
    core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership,
    biological_reaction_membership = expanded$reaction_membership,
    meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps,
    core_definition = paste(
      "complete Human-GEM GPR branch contained in the active Pando",
      "target-gene set for one condition-by-cell-type group"
    ),
    expansion_definition = paste(
      "one ordered pass: core subsystem, then KEGG/Reactome reaction",
      "equivalence, then master-Rhea reaction equivalence"
    ),
    analysis_group_unit =
      "condition_x_celltype_active_pando_metabolic_targets",
    feasibility_completion = "none_at_meta_module_stage"
  ))
  .rc_mm_write_tsv_gz(
    supported,
    file.path(outdir, "supported_metabolic_genes.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    core,
    file.path(outdir, "core_gene_reaction.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    expanded$reaction_membership,
    file.path(outdir, "meta_module_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    expanded$summary,
    file.path(outdir, "meta_module_summary.tsv.gz")
  )
  saveRDS(out, file.path(outdir, "condition_meta_modules.rds"))
  out
}
