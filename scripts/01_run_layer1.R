#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(RegCompassR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: 01_run_layer1.R <annotated_seurat.rds> <gpr_table.tsv> [output_dir]", call. = FALSE)
}
seurat_file <- args[[1]]
gpr_file <- args[[2]]
out_dir <- if (length(args) >= 3L) args[[3]] else "output/layer1"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object <- readRDS(seurat_file)
rc_validate_seurat_v4(object, rna_assay = "RNA", atac_assay = "ATAC", sample_col = "sample_id", celltype_col = "cell_type", condition_col = "condition")
inputs <- rc_extract_seurat_v4(object, rna_assay = "RNA", atac_assay = "ATAC", sample_col = "sample_id", celltype_col = "cell_type", condition_col = "condition")

pool_map <- rc_make_pools(inputs$meta, sample_col = "sample_id", celltype_col = "cell_type", condition_col = "condition", target_size = 80, min_pool_size = 30, min_group_size = 30, seed = 1)
rna_pb <- rc_pseudobulk_counts(inputs$rna_counts, pool_map, fun = "sum")
pool_meta <- rc_build_pool_metadata(pool_map, inputs$meta)
filtered <- rc_filter_empty_pools(rna_pb, pool_meta)
rna_logcpm <- rc_logcpm(filtered$counts)
pool_meta <- filtered$pool_meta
rna_detection <- rc_pool_detection(inputs$rna_counts, pool_map)
rna_detection <- rna_detection[, colnames(rna_logcpm), drop = FALSE]

gpr_table <- utils::read.delim(gpr_file, stringsAsFactors = FALSE)
layer1 <- rc_run_layer1_capacity(gpr_table = gpr_table, pool_expression = rna_logcpm, pool_detection = rna_detection, pool_meta = pool_meta, stratum_col = "cell_type", promiscuity_mode = "sqrt", and_method = "boltzmann", tau = 0.20)

rc_export_long_table(layer1$reaction_capacity_L1, file.path(out_dir, "layer1_capacity_long.tsv"), value_col = "C_rel")
utils::write.table(layer1$q95_diagnostics, file.path(out_dir, "q95_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(layer1$gpr_diagnostics, file.path(out_dir, "gpr_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rc_write_report_md(file.path(out_dir, "layer1_report.md"), q95_diagnostics = layer1$q95_diagnostics, gpr_diagnostics = layer1$gpr_diagnostics, confidence = layer1$reaction_confidence)
