.rc_condition_gene_regulatory_modifier_without_rsq_filter <-
  .rc_condition_gene_regulatory_modifier

.rc_condition_gene_regulatory_modifier <- function(
    significant_edges, object, unit_meta,
    condition_col = "condition", celltype_col = "cell_type",
    atac_assay = "ATAC", target_genes = NULL, min_scale = 0.05) {
  edges <- significant_edges
  if (is.data.frame(edges) && nrow(edges)) {
    rsq <- if ("rsq" %in% colnames(edges)) {
      suppressWarnings(as.numeric(edges$rsq))
    } else {
      rep(NA_real_, nrow(edges))
    }
    edges <- edges[is.finite(rsq), , drop = FALSE]
  }
  out <- .rc_condition_gene_regulatory_modifier_without_rsq_filter(
    significant_edges = edges,
    object = object,
    unit_meta = unit_meta,
    condition_col = condition_col,
    celltype_col = celltype_col,
    atac_assay = atac_assay,
    target_genes = target_genes,
    min_scale = min_scale
  )
  attr(out, "reliability_policy") <- paste(
    "only finite Pando R-squared values are trusted; targets without finite",
    "R-squared receive regulatory reliability zero"
  )
  out
}

.rc_feasibility_completion_metadata <- function(model_mode) {
  if (identical(model_mode, "meta_module_gem")) {
    return(list(
      feasibility_completion =
        "local_unconstrained_fastcore_then_global_union_medium_specific_fastcore",
      feasibility_completion_stages = list(
        local = paste(
          "condition x cell-type biological meta-modules are completed against",
          "an unconstrained shared FASTCC parent"
        ),
        global = paste(
          "the union model is rebuilt and add-only FASTCORE-completed separately",
          "for each shared medium scenario before scoring"
        )
      )
    ))
  }
  list(
    feasibility_completion = "not_applicable_full_gem",
    feasibility_completion_stages = list(
      local = "local meta-modules are not used for full-GEM scoring",
      global = "the complete GEM is constrained by each shared medium"
    )
  )
}
