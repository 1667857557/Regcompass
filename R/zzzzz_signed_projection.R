# Signed projection policy loaded after the metadata-enrichment wrapper.
# Discordant or mixed shared-TF relations remain in diagnostics but do not join
# genes into one biological component. Direct TF -> target relations remain
# valid component edges regardless of activation or repression sign.

.rc_signed_metadata_project_metabolic_grn <- rc_project_metabolic_grn

rc_project_metabolic_grn <- function(tf_peak_gene, metabolic_genes,
                                     top_k = 5L, min_shared_tfs = 1L,
                                     min_tf_jaccard = 0,
                                     max_targets_per_tf = 200L,
                                     include_direct_metabolic_tf = TRUE) {
  answer <- .rc_signed_metadata_project_metabolic_grn(
    tf_peak_gene = tf_peak_gene,
    metabolic_genes = metabolic_genes,
    top_k = top_k,
    min_shared_tfs = min_shared_tfs,
    min_tf_jaccard = min_tf_jaccard,
    max_targets_per_tf = max_targets_per_tf,
    include_direct_metabolic_tf = include_direct_metabolic_tf
  )
  edges <- answer$edges
  nodes <- answer$nodes
  if (!nrow(edges) || !nrow(nodes)) {
    edges$used_for_component <- logical(nrow(edges))
    answer$edges <- edges
    return(answer)
  }

  relation <- as.character(edges$regulatory_relation)
  edges$used_for_component <- edges$direct_regulatory %in% TRUE |
    (!edges$direct_regulatory %in% TRUE & relation == "concordant")
  edges$used_for_component[is.na(edges$used_for_component)] <- FALSE

  for (sample in unique(as.character(nodes$sample_id))) {
    node_index <- as.character(nodes$sample_id) == sample
    edge_index <- as.character(edges$sample_id) == sample
    component_edges <- edges[
      edge_index & edges$used_for_component,
      , drop = FALSE
    ]
    components <- .rc_mm_components(
      as.character(nodes$gene[node_index]),
      component_edges
    )
    module_ids <- paste0(
      sample, "::GRN", sprintf("%04d", components$component)
    )
    nodes$module_id[node_index] <- module_ids[
      match(as.character(nodes$gene[node_index]), components$gene)
    ]
    if (any(edge_index)) {
      edges$module_id[edge_index] <- nodes$module_id[
        match(
          as.character(edges$gene_a[edge_index]),
          as.character(nodes$gene)
        )
      ]
    }
  }

  answer$nodes <- nodes
  answer$edges <- edges
  answer$component_policy <- paste(
    "direct regulatory edges plus concordant signed shared-TF edges;",
    "discordant and mixed shared-TF edges are diagnostic only"
  )
  answer
}
