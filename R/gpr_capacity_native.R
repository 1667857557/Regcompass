.rc_compile_gpr_indices <- function(gpr_list, gene_ids) {
  reaction_group_offset <- integer(length(gpr_list) + 1L)
  group_gene_offset <- 0L
  gene_index <- integer()
  group_count <- 0L

  for (reaction in seq_along(gpr_list)) {
    groups <- gpr_list[[reaction]]
    for (group in groups) {
      genes <- unique(tolower(group))
      index <- match(genes, gene_ids)
      index[is.na(index)] <- 0L
      gene_index <- c(gene_index, as.integer(index))
      group_gene_offset <- c(
        group_gene_offset, as.integer(length(gene_index))
      )
      group_count <- group_count + 1L
    }
    reaction_group_offset[[reaction + 1L]] <- group_count
  }

  list(
    reaction_group_offset = as.integer(reaction_group_offset),
    group_gene_offset = as.integer(group_gene_offset),
    gene_index = as.integer(gene_index)
  )
}
