.rc_summarize_supported_metabolic_genes <- function(
    grn_result, metabolic_genes) {
  if (!is.list(grn_result) ||
      !is.data.frame(grn_result$tf_peak_gene_significant)) {
    stop("`grn_result` is not a valid single-cell GRN result.", call. = FALSE)
  }
  group_cols <- as.character(grn_result$group_cols)
  required <- c("group_id", group_cols, "tf", "target", "region")
  significant <- grn_result$tf_peak_gene_significant
  missing <- setdiff(required, colnames(significant))
  if (length(missing)) {
    stop(
      "Active GRN edge table is missing columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  metabolic_genes <- unique(toupper(.rc_mm_trim_unique(metabolic_genes)))
  significant$target <- toupper(trimws(as.character(significant$target)))
  significant$tf <- toupper(trimws(as.character(significant$tf)))
  significant$region <- trimws(as.character(significant$region))
  significant <- significant[
    significant$target %in% metabolic_genes,
    , drop = FALSE
  ]
  if (!nrow(significant)) {
    stop(
      "No GEM metabolic target genes have active TF-peak-gene evidence.",
      call. = FALSE
    )
  }

  key <- paste(significant$group_id, significant$target, sep = "\001")
  rows <- split(seq_len(nrow(significant)), key)
  summary_rows <- lapply(rows, function(index) {
    one <- significant[index, , drop = FALSE]
    group_id <- as.character(one$group_id[[1L]])
    target <- as.character(one$target[[1L]])
    numeric_column <- function(name) {
      if (name %in% colnames(one)) {
        suppressWarnings(as.numeric(one[[name]]))
      } else {
        rep(NA_real_, nrow(one))
      }
    }
    estimate <- numeric_column("estimate")
    effective <- if ("effective_estimate" %in% colnames(one)) {
      numeric_column("effective_estimate")
    } else {
      estimate
    }
    padj <- numeric_column("padj")
    rsq <- if ("cv_rsq" %in% colnames(one)) {
      numeric_column("cv_rsq")
    } else {
      numeric_column("rsq")
    }
    selection <- numeric_column("selection_frequency")
    sign_stability <- numeric_column("sign_stability")
    finite_min <- function(value) {
      value <- value[is.finite(value)]
      if (length(value)) min(value) else NA_real_
    }
    finite_max <- function(value) {
      value <- value[is.finite(value)]
      if (length(value)) max(value) else NA_real_
    }
    evidence <- if ("evidence_type" %in% colnames(one)) {
      paste(sort(unique(as.character(one$evidence_type))), collapse = ";")
    } else {
      "legacy_significant_pando_tf_peak_gene_target"
    }
    group_values <- one[1L, group_cols, drop = FALSE]
    data.frame(
      group_id = group_id,
      group_values,
      sample_id = group_id,
      module_id = paste0(group_id, "::SUPPORTED_METABOLIC_GENES"),
      gene = target,
      n_significant_edges = nrow(one),
      n_regulating_tfs = length(unique(one$tf[nzchar(one$tf)])),
      n_regulatory_regions = length(unique(one$region[nzchar(one$region)])),
      min_padj = finite_min(padj),
      max_abs_estimate = finite_max(abs(estimate)),
      max_abs_effective_estimate = finite_max(abs(effective)),
      max_model_rsq = finite_max(rsq),
      min_selection_frequency = finite_min(selection),
      min_sign_stability = finite_min(sign_stability),
      n_positive_edges = sum(is.finite(effective) & effective > 0),
      n_negative_edges = sum(is.finite(effective) & effective < 0),
      evidence_definition = evidence,
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
        "Active GRN-supported metabolic genes did not completely satisfy",
        "any GEM GPR branch."
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

  multitask <- identical(
    as.character(grn_result$grn_mode %||% ""),
    "multitask_shared_backbone"
  )
  out <- c(grn_result, list(
    supported_metabolic_genes = supported,
    core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership,
    biological_reaction_membership = expanded$reaction_membership,
    meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps,
    core_definition = if (multitask) {
      paste(
        "complete GEM GPR branch contained in the stability-selected",
        "condition sub-GRN target-gene set for one condition-by-cell-type group"
      )
    } else {
      paste(
        "complete GEM GPR branch contained in the significant Pando",
        "target-gene set for one condition-by-cell-type group"
      )
    },
    expansion_definition = paste(
      "one ordered pass: core subsystem, then KEGG/Reactome reaction",
      "equivalence, then master-Rhea reaction equivalence"
    ),
    analysis_group_unit = if (multitask) {
      "condition_x_celltype_stability_selected_multitask_metabolic_targets"
    } else {
      "condition_x_celltype_significant_pando_metabolic_targets"
    },
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
