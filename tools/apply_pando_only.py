from pathlib import Path
import re


def function_end(text, start):
    brace = text.index("{", start)
    depth = 0
    quote = None
    escaped = False
    for i in range(brace, len(text)):
        ch = text[i]
        if quote is not None:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == quote:
                quote = None
            continue
        if ch in ('"', "'"):
            quote = ch
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i + 1
    raise RuntimeError("Unbalanced R function")


global_path = Path("R/global_workflow.R")
global_text = global_path.read_text()
stratum_start = global_text.index(".rc_run_regcompass_stratum <- function")
stratum_end = function_end(global_text, stratum_start)
stratum_replacement = r'''.rc_metacell_logcpm <- function(counts, scale_factor = 1e6) {
  counts <- methods::as(counts, "dgCMatrix")
  library_size <- Matrix::colSums(counts)
  if (any(!is.finite(library_size)) || any(library_size <= 0)) {
    stop("Every metacell must have a positive finite RNA library size.", call. = FALSE)
  }
  scaled <- counts %*% Matrix::Diagonal(x = scale_factor / library_size)
  log1p(scaled)
}

.rc_run_regcompass_stratum <- function(object, group_id, group_cols, gem, outdir,
                                        pfm, genome, fragment_files = NULL,
                                        sample_col = "sample_id",
                                        condition_col = "condition",
                                        celltype_col = "cell_type",
                                        rna_assay = "RNA",
                                        atac_assay = "ATAC",
                                        metacell_args = list(),
                                        layer1_args = list(),
                                        pando_args = list()) {
  retired_link_args <- intersect(
    names(layer1_args),
    c(
      "linkpeaks_args", "min_metacells_for_linkpeaks",
      "force_metacell_relink", "allow_supplied_links",
      "peak_gene_links", "recompute_peak_gene_links"
    )
  )
  if (length(retired_link_args)) {
    stop(
      paste0(
        "Integrated RegCompass no longer runs a separate LinkPeaks step; ",
        "Pando performs peak-gene modeling internally. Remove: ",
        paste(retired_link_args, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  meta <- object@meta.data
  ids <- rc_make_stratum_id(meta, group_cols)
  cells <- rownames(meta)[ids == group_id]
  if (!length(cells)) {
    stop("No cells found for stratum: ", group_id, call. = FALSE)
  }
  one <- subset(object, cells = cells)
  stratum_dir <- file.path(
    outdir,
    gsub("[^A-Za-z0-9_.-]+", "_", group_id)
  )
  dir.create(stratum_dir, recursive = TRUE, showWarnings = FALSE)
  capacity_params <- list(
    promiscuity_mode = layer1_args$promiscuity_mode %||% "sqrt",
    and_method = layer1_args$and_method %||% "boltzmann",
    tau = layer1_args$tau %||% 0.20,
    or_method = "sum_sqrtK"
  )
  metacell_defaults <- list(
    object = one,
    outdir = file.path(stratum_dir, "01_metacells"),
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    fragment_files = fragment_files,
    save_metacell_object = TRUE,
    save_counts = TRUE,
    save_fragments = TRUE,
    require_fragment_aggregation = TRUE,
    fragment_aggregation_backend = "regcompass",
    BPPARAM = FALSE,
    on_stratum_error = "stop"
  )
  reserved <- intersect(names(metacell_args), names(metacell_defaults))
  reserved <- setdiff(
    reserved,
    c(
      "rna_reduction", "atac_reduction", "rna_dims", "atac_dims",
      "gamma", "seed", "min_cells_per_stratum",
      "min_metacell_size", "min_metacells_per_stratum",
      "adaptive_gamma", "label_col", "bgzip_path", "tabix_path",
      "fragment_nb_cl", "overwrite"
    )
  )
  if (length(reserved)) {
    stop(
      "`metacell_args` cannot override workflow fields: ",
      paste(reserved, collapse = ", "),
      call. = FALSE
    )
  }
  metacell_defaults[names(metacell_args)] <- NULL
  metacells <- do.call(
    rc_make_supercell2_metacells,
    c(metacell_defaults, metacell_args)
  )
  minimum_metacells <- as.integer(pando_args$min_metacells %||% 20L)
  if (nrow(metacells$metacell_meta) < minimum_metacells) {
    stop(
      "Stratum `", group_id, "` produced fewer than ",
      minimum_metacells, " metacells.",
      call. = FALSE
    )
  }
  metacell_object <- rc_load_or_merge_metacell_objects(
    metacells$metacell_objects,
    fragment_manifest = metacells$fragment_manifest,
    metacell_meta = metacells$metacell_meta,
    fragment_files = metacells$fragment_files,
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    require_complete_fragments = TRUE
  )
  pando_args$BPPARAM <- NULL
  pando_args$save_sample_metacell_objects <- NULL
  pando_infer_args <- pando_args$pando_infer_args %||% list()
  pando_infer_args$parallel <- FALSE
  pando_args$pando_infer_args <- pando_infer_args
  pando_args$group_cols <- NULL
  pando_args$sample_col <- NULL
  pando_args$condition_col <- NULL
  pando_args$celltype_col <- NULL
  pando_outdir <- file.path(stratum_dir, "02_pando_meta_modules")
  pando_defaults <- list(
    metacell_object = metacell_object,
    gem = gem,
    outdir = pando_outdir,
    pfm = pfm,
    genome = genome,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    group_cols = group_cols,
    single_cell_genes = rownames(.rc_get_assay_counts(one, rna_assay)),
    rna_assay = rna_assay,
    atac_assay = atac_assay,
    min_metacells = minimum_metacells,
    save_sample_metacell_objects = TRUE,
    BPPARAM = FALSE,
    on_sample_error = "stop"
  )
  pando_defaults[names(pando_args)] <- NULL
  meta_modules <- do.call(
    rc_run_pando_meta_modules,
    c(pando_defaults, pando_args)
  )
  pando_object_file <- file.path(
    pando_outdir,
    "sample_metacell_objects",
    paste0(gsub("[^A-Za-z0-9_.-]+", "_", group_id), ".rds")
  )
  if (!file.exists(pando_object_file)) {
    stop(
      "Pando did not save its normalized metacell object: ",
      pando_object_file,
      call. = FALSE
    )
  }
  pando_object <- readRDS(pando_object_file)
  pando_confidence <- .rc_pando_reaction_confidence(
    meta_modules,
    pando_object,
    gem,
    atac_assay = atac_assay
  )
  saveRDS(
    pando_confidence$gene_confidence,
    file.path(pando_outdir, "pando_gene_confidence.rds")
  )
  saveRDS(
    pando_confidence$reaction_confidence_matrix,
    file.path(pando_outdir, "pando_reaction_confidence_matrix.rds")
  )
  .rc_mm_write_tsv_gz(
    pando_confidence$gene_confidence_diagnostics,
    file.path(pando_outdir, "pando_gene_confidence_diagnostics.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    pando_confidence$reaction_confidence,
    file.path(pando_outdir, "pando_reaction_confidence.tsv.gz")
  )
  meta_modules$gene_confidence <- pando_confidence$gene_confidence
  meta_modules$gene_confidence_diagnostics <-
    pando_confidence$gene_confidence_diagnostics
  meta_modules$reaction_confidence <-
    pando_confidence$reaction_confidence
  meta_modules$reaction_confidence_matrix <-
    pando_confidence$reaction_confidence_matrix
  meta_modules$confidence_source <- pando_confidence$confidence_source

  gpr_genes <- toupper(rc_metabolic_gpr_genes(gem$gpr_table))
  rna_counts <- metacells$rna_counts[
    toupper(rownames(metacells$rna_counts)) %in% gpr_genes,
    ,
    drop = FALSE
  ]
  if (!nrow(rna_counts)) {
    stop(
      "No Human-GEM GPR genes were retained in metacell RNA counts.",
      call. = FALSE
    )
  }
  rna_logcpm <- .rc_metacell_logcpm(rna_counts)
  unit_meta <- metacells$metacell_meta
  id_col <- if ("metacell_id" %in% colnames(unit_meta)) {
    "metacell_id"
  } else if ("pool_id" %in% colnames(unit_meta)) {
    "pool_id"
  } else {
    NULL
  }
  if (is.null(id_col)) {
    stop("Metacell metadata lacks metacell_id/pool_id.", call. = FALSE)
  }
  if (!"pool_id" %in% colnames(unit_meta)) {
    unit_meta$pool_id <- as.character(unit_meta[[id_col]])
  }
  if (!"unit_id" %in% colnames(unit_meta)) {
    unit_meta$unit_id <- unit_meta$pool_id
  }
  unit_meta <- unit_meta[
    match(colnames(rna_logcpm), as.character(unit_meta$pool_id)),
    ,
    drop = FALSE
  ]
  if (anyNA(unit_meta$pool_id)) {
    stop("Metacell metadata are incomplete after Pando.", call. = FALSE)
  }
  reaction_confidence <- as.matrix(
    pando_confidence$reaction_confidence_matrix
  )
  missing_confidence_units <- setdiff(
    colnames(rna_logcpm),
    colnames(reaction_confidence)
  )
  if (length(missing_confidence_units)) {
    stop(
      "Pando reaction confidence is missing metacells: ",
      paste(utils::head(missing_confidence_units, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  reaction_confidence <- reaction_confidence[
    ,
    colnames(rna_logcpm),
    drop = FALSE
  ]
  layer1 <- list(
    schema_version = "regcompass_stratum_evidence_v2",
    rna_metacell_logcpm = rna_logcpm,
    reaction_confidence = reaction_confidence,
    reaction_confidence_source =
      "pando_internal_peak_gene_accessibility",
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    strict_group_id = group_id
  )
  artifact <- list(
    schema_version = "regcompass_stratum_v2",
    group_id = group_id,
    group_cols = group_cols,
    capacity_params = capacity_params,
    layer1 = layer1,
    grn_meta_modules = meta_modules,
    metacell_meta = unit_meta,
    metacell_dir = stratum_dir
  )
  artifact_file <- file.path(stratum_dir, "stratum_result.rds")
  saveRDS(artifact, artifact_file)
  list(
    group_id = group_id,
    status = "ok",
    artifact_file = artifact_file,
    n_cells = length(cells),
    n_metacells = nrow(unit_meta),
    error_class = NA_character_,
    error_message = NA_character_
  )
}'''
global_text = (
    global_text[:stratum_start]
    + stratum_replacement
    + global_text[stratum_end:]
)

merge_start = global_text.index(".rc_merge_stratum_layer1 <- function")
merge_end = function_end(global_text, merge_start)
merge_text = global_text[merge_start:merge_end]
pattern = re.compile(
    r"  confidence_list <- lapply\(artifacts, function\(x\) \{.*?"
    r"  reaction_confidence <- \.rc_cbind_matrix_union\(confidence_list\)",
    re.S,
)
replacement = '''  confidence_list <- lapply(artifacts, function(x) {
    value <- x$layer1$reaction_confidence
    if (is.null(value)) {
      stop(
        "Every upstream artifact must contain Pando-derived reaction confidence.",
        call. = FALSE
      )
    }
    as.matrix(value)
  })
  reaction_confidence <- .rc_cbind_matrix_union(confidence_list)'''
merge_text, count = pattern.subn(replacement, merge_text, count=1)
if count != 1:
    raise RuntimeError("Global confidence merge block was not replaced")
merge_text = merge_text.replace(
    'capacity_calibration_scope = "all_metacells_global_gene_score_and_reaction_q95",',
    'capacity_calibration_scope = "all_metacells_global_gene_score_and_reaction_q95",\n'
    '    reaction_confidence_source = "pando_internal_peak_gene_accessibility",',
)
global_text = global_text[:merge_start] + merge_text + global_text[merge_end:]
global_path.write_text(global_text)


test_path = Path("tests/testthat/test_global_metacell_workflow.R")
test_text = test_path.read_text()
calibration_start = test_text.index(
    'test_that("global Layer 1 recomputes gene scores and reaction Q95 across all metacells"'
)
calibration_end = function_end(test_text, calibration_start)
if calibration_end < len(test_text) and test_text[calibration_end] == ")":
    calibration_end += 1
while calibration_end < len(test_text) and test_text[calibration_end] in "\r\n":
    calibration_end += 1
calibration_test = r'''test_that("global Layer 1 recomputes gene scores and uses Pando confidence", {
  artifacts <- list(
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      layer1 = list(
        rna_metacell_logcpm = matrix(
          1, nrow = 1, dimnames = list("G1", "u1")
        ),
        reaction_confidence = matrix(
          0.8, nrow = 1, dimnames = list("R1", "u1")
        ),
        unit_meta = data.frame(
          pool_id = "u1", unit_id = "u1", sample_id = "S1",
          condition = "A", cell_type = "T", stringsAsFactors = FALSE
        )
      )
    ),
    list(
      capacity_params = list(
        promiscuity_mode = "sqrt", and_method = "boltzmann",
        tau = 0.20, or_method = "sum_sqrtK"
      ),
      layer1 = list(
        rna_metacell_logcpm = matrix(
          3, nrow = 1, dimnames = list("G1", "u2")
        ),
        reaction_confidence = matrix(
          0.4, nrow = 1, dimnames = list("R1", "u2")
        ),
        unit_meta = data.frame(
          pool_id = "u2", unit_id = "u2", sample_id = "S2",
          condition = "B", cell_type = "T", stringsAsFactors = FALSE
        )
      )
    )
  )
  gem <- list(gpr_table = data.frame(
    reaction_id = "R1", and_group_id = 1, gene = "G1"
  ))
  out <- .rc_merge_stratum_layer1(
    artifacts, gem, single_cell_genes = "G1",
    sample_col = "sample_id", condition_col = "condition",
    celltype_col = "cell_type"
  )
  expect_equal(
    out$capacity_calibration_scope,
    "all_metacells_global_gene_score_and_reaction_q95"
  )
  expect_equal(
    out$reaction_confidence_source,
    "pando_internal_peak_gene_accessibility"
  )
  expect_equal(unique(as.character(out$q95_diagnostics$stratum)), "global")
  expect_identical(colnames(out$C_rel), c("u1", "u2"))
  expect_equal(out$reaction_confidence["R1", ], c(u1 = 0.8, u2 = 0.4))
  expect_lt(out$global_gene_score["G1", "u1"], out$global_gene_score["G1", "u2"])
  expect_lt(out$C_raw["R1", "u1"], out$C_raw["R1", "u2"])
})

'''
test_text = (
    test_text[:calibration_start]
    + calibration_test
    + test_text[calibration_end:]
)
insert_anchor = (
    'test_that("global meta-module union preserves source tables and creates one canonical module", {'
)
no_link_test = r'''test_that("integrated stratum worker uses Pando without rerunning LinkPeaks", {
  body_text <- paste(
    deparse(body(.rc_run_regcompass_stratum)),
    collapse = "\n"
  )
  expect_match(body_text, "rc_run_pando_meta_modules", fixed = TRUE)
  expect_match(
    body_text,
    "pando_internal_peak_gene_accessibility",
    fixed = TRUE
  )
  expect_false(grepl("rc_run_layer1_from_metacells", body_text, fixed = TRUE))
  expect_false(grepl("Signac::LinkPeaks", body_text, fixed = TRUE))
})

'''
if no_link_test not in test_text:
    test_text = test_text.replace(
        insert_anchor,
        no_link_test + insert_anchor,
        1,
    )
test_path.write_text(test_text)


readme_path = Path("README.md")
readme = readme_path.read_text()
readme = readme.replace(
    "     metacell LinkPeaks\n"
    "     Layer 1 RNA-GPR capacity and ATAC confidence\n"
    "     Pando GRN",
    "     Pando internal peak-gene and TF-gene modeling\n"
    "     Pando-derived metacell ATAC confidence\n"
    "     Pando GRN",
)
readme = readme.replace(
    "The shared GEM is the structural reference for every metacell. "
    "Biological differences enter the LP objective through metacell-specific "
    "penalties, not through sample-specific stoichiometric models. The final "
    "expression-capacity normalization recomputes gene scores from the combined "
    "GPR-gene logCPM matrix and then applies one global reaction-wise Q95; "
    "stratum-local expression capacities are not used for cross-sample scoring.",
    "The shared GEM is the structural reference for every metacell. "
    "Biological differences enter the LP objective through metacell-specific "
    "penalties, not through sample-specific stoichiometric models. The final "
    "expression-capacity normalization recomputes gene scores from the combined "
    "GPR-gene logCPM matrix and then applies one global reaction-wise Q95. "
    "Peak-gene evidence is computed once inside Pando; RegCompass derives "
    "metacell ATAC confidence from significant Pando regions and does not run "
    "a separate LinkPeaks pass.",
)
readme = re.sub(
    r"  layer1_args = list\(\n"
    r"    min_metacells_for_linkpeaks = 10,\n"
    r"    bootstrap = FALSE\n"
    r"  \),",
    '  layer1_args = list(\n'
    '    promiscuity_mode = "sqrt",\n'
    '    and_method = "boltzmann",\n'
    '    tau = 0.20\n'
    '  ),',
    readme,
)
readme_path.write_text(readme)


news_path = Path("NEWS.md")
news = news_path.read_text()
note = (
    "- Removed the duplicate Signac LinkPeaks pass from the integrated workflow; "
    "Pando is now the sole peak-gene model, and significant Pando regions plus "
    "metacell accessibility define ATAC confidence.\n"
)
if note not in news:
    first_header_end = news.find("\n", news.find("# "))
    news = news[: first_header_end + 1] + "\n" + note + news[first_header_end + 1 :]
news_path.write_text(news)


Path(".github/workflows/apply-pando-only-peak-gene.yml").unlink(missing_ok=True)
Path("tools/apply_pando_only.py").unlink(missing_ok=True)
