from pathlib import Path
import re


def function_span(text, name):
    pattern = re.compile(r'(?m)^' + re.escape(name) + r'\s*<-\s*function\b')
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f"{name}: expected one definition, found {len(matches)}")
    start = matches[0].start()
    opening = text.find("{", matches[0].end())
    depth = 0
    quote = None
    escaped = False
    comment = False
    for i in range(opening, len(text)):
        ch = text[i]
        if comment:
            if ch == "\n":
                comment = False
            continue
        if quote:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == quote:
                quote = None
            continue
        if ch == "#":
            comment = True
        elif ch in ("'", '"'):
            quote = ch
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return start, i + 1
    raise RuntimeError(f"{name}: closing brace not found")


def replace_function(path, name, replacement):
    p = Path(path)
    text = p.read_text()
    start, end = function_span(text, name)
    p.write_text(text[:start] + replacement.rstrip() + text[end:])


def replace_once(path, old, new, label):
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    p.write_text(text.replace(old, new, 1))


stage_filter = r'''.rc_filter_stage1_groups_by_min_cells <- function(
    object, condition_col, celltype_col, cell_type, min_cells) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from Seurat.", call. = FALSE)
  }
  .rc_validate_celltype_metadata(object@meta.data, celltype_col)
  observed_type <- trimws(as.character(object@meta.data[[celltype_col]]))
  available_type <- unique(observed_type)
  requested_type <- if (is.null(cell_type)) {
    available_type
  } else {
    unique(trimws(as.character(cell_type)))
  }
  missing_type <- setdiff(requested_type, available_type)
  if (length(missing_type)) {
    stop("Requested cell types were not found: ",
         paste(missing_type, collapse = ", "), call. = FALSE)
  }
  selected <- observed_type %in% requested_type
  if (!any(selected)) {
    stop("No cells remain after applying `cell_type`.", call. = FALSE)
  }

  has_condition <- is.character(condition_col) && length(condition_col) == 1L &&
    !is.na(condition_col) && nzchar(condition_col) &&
    condition_col %in% colnames(object@meta.data)
  if (has_condition) {
    observed_condition <- trimws(as.character(object@meta.data[[condition_col]]))
    invalid_condition <- selected &
      (is.na(object@meta.data[[condition_col]]) | !nzchar(observed_condition))
    if (any(invalid_condition)) {
      stop("Condition metadata contain missing or empty values in requested cell types.",
           call. = FALSE)
    }
    condition_levels <- unique(observed_condition[selected])
  } else {
    observed_condition <- rep(NA_character_, nrow(object@meta.data))
    condition_levels <- character()
  }

  condition_design <- length(condition_levels) >= 2L
  if (condition_design) {
    count_matrix <- table(
      factor(observed_type[selected], levels = requested_type),
      factor(observed_condition[selected], levels = condition_levels)
    )
    retained_stratum <- count_matrix >= as.integer(min_cells)
    retained_condition_count <- base::rowSums(retained_stratum)
    retained_type <- rownames(count_matrix)[retained_condition_count >= 1L]
    condition_type <- rownames(count_matrix)[retained_condition_count >= 2L]
    standard_type <- rownames(count_matrix)[retained_condition_count == 1L]

    diagnostics <- do.call(rbind, lapply(requested_type, function(type) {
      stratum_ok <- as.logical(retained_stratum[type, condition_levels])
      type_mode <- if (type %in% condition_type) {
        "condition_grn"
      } else if (type %in% standard_type) {
        "standard_pando"
      } else {
        "excluded"
      }
      data.frame(
        cell_type = type,
        condition = condition_levels,
        n_cells = as.integer(count_matrix[type, condition_levels]),
        retained_stratum = stratum_ok,
        retained_cell_type = type %in% retained_type,
        cell_type_analysis_mode = type_mode,
        retained = stratum_ok,
        fit_status = ifelse(
          !stratum_ok,
          "excluded_below_min_cells",
          ifelse(type %in% condition_type,
                 "eligible_condition_pando",
                 "eligible_standard_pando_single_condition")
        ),
        n_retained_conditions = as.integer(retained_condition_count[[type]]),
        threshold = as.integer(min_cells),
        threshold_scope = "condition_x_cell_type",
        stringsAsFactors = FALSE
      )
    }))

    dropped <- diagnostics[!diagnostics$retained_stratum, , drop = FALSE]
    if (nrow(dropped)) {
      message(
        "Stage 1 excluded condition x cell-type strata below min_cells=300: ",
        paste0(dropped$cell_type, "{", dropped$condition, "}=",
               dropped$n_cells, collapse = "; ")
      )
    }
    if (!length(retained_type)) {
      stop("No condition x cell-type stratum reaches min_cells=300.",
           call. = FALSE)
    }
    if (length(standard_type)) {
      message(
        "Using standard Pando for cell types with one retained condition: ",
        paste(standard_type, collapse = ", ")
      )
    }

    keep_cells <- rep(FALSE, nrow(object@meta.data))
    selected_rows <- which(selected)
    type_index <- match(observed_type[selected_rows], rownames(retained_stratum))
    condition_index <- match(
      observed_condition[selected_rows], colnames(retained_stratum)
    )
    valid <- !is.na(type_index) & !is.na(condition_index)
    selected_keep <- rep(FALSE, length(selected_rows))
    selected_keep[valid] <- retained_stratum[cbind(
      type_index[valid], condition_index[valid]
    )]
    keep_cells[selected_rows] <- selected_keep
    analysis_mode <- if (length(condition_type) && length(standard_type)) {
      "mixed_pando"
    } else if (length(condition_type)) {
      "condition_grn"
    } else {
      "standard_pando"
    }
  } else {
    counts <- table(factor(observed_type[selected], levels = requested_type))
    retained_type <- names(counts)[as.integer(counts) >= min_cells]
    condition_type <- character()
    standard_type <- retained_type
    diagnostics <- data.frame(
      cell_type = names(counts),
      condition = if (length(condition_levels) == 1L) {
        condition_levels[[1L]]
      } else {
        NA_character_
      },
      n_cells = as.integer(counts),
      retained_stratum = as.integer(counts) >= min_cells,
      retained_cell_type = names(counts) %in% retained_type,
      cell_type_analysis_mode = ifelse(
        names(counts) %in% retained_type, "standard_pando", "excluded"
      ),
      retained = names(counts) %in% retained_type,
      fit_status = ifelse(
        names(counts) %in% retained_type,
        "eligible_standard_pando",
        "excluded_below_min_cells"
      ),
      n_retained_conditions = if (length(condition_levels) == 1L) {
        as.integer(names(counts) %in% retained_type)
      } else {
        NA_integer_
      },
      threshold = as.integer(min_cells),
      threshold_scope = "cell_type",
      stringsAsFactors = FALSE
    )
    if (!length(retained_type)) {
      stop("No requested cell type reaches min_cells=300.", call. = FALSE)
    }
    keep_cells <- selected & observed_type %in% retained_type
    analysis_mode <- "standard_pando"
  }

  retained_cells <- rownames(object@meta.data)[keep_cells]
  filtered <- subset(object, cells = retained_cells)
  filtered@misc$regcompass_stage1_group_filter <- diagnostics
  list(
    object = filtered,
    retained_cell_types = retained_type,
    condition_pando_cell_types = condition_type,
    standard_pando_cell_types = standard_type,
    diagnostics = diagnostics,
    analysis_mode = analysis_mode,
    condition_levels = if (condition_design) {
      unique(observed_condition[keep_cells])
    } else {
      condition_levels
    }
  )
}'''

stage_build = r'''.rc_build_stage_analysis_cell_set <- function(
    object, condition_col, celltype_col, cell_type, pando_args) {
  threshold <- .rc_resolve_stage1_min_cells_contract(pando_args)
  groups <- .rc_filter_stage1_groups_by_min_cells(
    object = object,
    condition_col = condition_col,
    celltype_col = celltype_col,
    cell_type = cell_type,
    min_cells = threshold$min_cells
  )
  analysis_types <- unique(as.character(groups$retained_cell_types))
  observed_type <- trimws(as.character(groups$object@meta.data[[celltype_col]]))
  analysis_cells <- rownames(groups$object@meta.data)[
    observed_type %in% analysis_types
  ]
  if (!length(analysis_cells)) {
    stop("No cells remain for Stage 1 Pando analysis.", call. = FALSE)
  }
  list(
    object = subset(groups$object, cells = analysis_cells),
    retained_cells = analysis_cells,
    retained_cell_types = analysis_types,
    condition_pando_cell_types =
      unique(as.character(groups$condition_pando_cell_types)),
    standard_pando_cell_types =
      unique(as.character(groups$standard_pando_cell_types)),
    diagnostics = groups$diagnostics,
    analysis_mode = groups$analysis_mode,
    condition_levels = groups$condition_levels,
    min_cells = threshold$min_cells,
    pando_args = threshold$pando_args
  )
}'''

replace_function("R/stage1_input_contract.R",
                 ".rc_filter_stage1_groups_by_min_cells", stage_filter)
replace_function("R/stage1_input_contract.R",
                 ".rc_build_stage_analysis_cell_set", stage_build)

mixed_file = r'''# Cell-type-specific routing between common-dictionary and standard Pando.

.rc_bind_pando_field <- function(results, field) {
  values <- lapply(results, function(x) x[[field]])
  values <- values[vapply(values, is.data.frame, logical(1))]
  values <- values[vapply(values, nrow, integer(1)) > 0L]
  if (!length(values)) return(data.frame())
  .rc_bind_frames_fill(values)
}

.rc_standard_pando_runtime_args <- function(args) {
  if (!is.list(args)) return(list())
  method <- tryCatch(
    getS3method("infer_grn", "GRNData", envir = asNamespace("Pando")),
    error = function(error) NULL
  )
  if (!is.function(method)) return(list())
  allowed <- setdiff(names(formals(method)), c("object", "..."))
  args[intersect(names(args), allowed)]
}

.rc_merge_pando_results <- function(
    condition_result = NULL, standard_results = list(),
    condition_types = character(), standard_types = character(),
    condition_col, celltype_col, outdir) {
  results <- c(
    if (is.null(condition_result)) list() else list(condition_result),
    standard_results
  )
  if (!length(results)) stop("No Pando result was produced.", call. = FALSE)
  standard_objects <- list()
  for (value in standard_results) {
    for (name in names(value$standard_pando_objects)) {
      if (name %in% names(standard_objects)) {
        stop("Duplicated standard-Pando cell type: ", name, call. = FALSE)
      }
      standard_objects[[name]] <- value$standard_pando_objects[[name]]
    }
  }
  condition_fits <- if (is.null(condition_result)) list() else
    condition_result$condition_grn_fits
  paired_meta <- .rc_bind_pando_field(results, "paired_cell_metadata")
  if (nrow(paired_meta)) {
    paired_meta <- paired_meta[!duplicated(paired_meta$cell_id), , drop = FALSE]
  }
  routing <- rbind(
    if (length(condition_types)) data.frame(
      cell_type = condition_types, analysis_mode = "condition_grn",
      stringsAsFactors = FALSE
    ),
    if (length(standard_types)) data.frame(
      cell_type = standard_types, analysis_mode = "standard_pando",
      stringsAsFactors = FALSE
    )
  )
  rownames(routing) <- NULL
  mode <- if (length(condition_types) && length(standard_types)) {
    "mixed_pando"
  } else if (length(condition_types)) {
    "condition_grn"
  } else {
    "standard_pando"
  }
  answer <- list(
    schema_version = "regcompass_celltype_routed_pando_v1",
    analysis_mode = mode,
    cell_type_analysis_mode = routing,
    condition_coefficients_calculated = length(condition_fits) > 0L,
    pando_fit_schema = if (length(condition_fits)) {
      .RC_PANDO_CONDITION_GRN_FIT_SCHEMA
    } else NA_character_,
    pando_installed_version = as.character(utils::packageVersion("Pando")),
    pando_grn_data = if (is.null(condition_result)) NULL else
      condition_result$pando_grn_data,
    standard_pando_objects = standard_objects,
    condition_grn_fits = condition_fits,
    target_metabolic_genes = unique(unlist(
      lapply(results, `[[`, "target_metabolic_genes"), use.names = FALSE
    )),
    condition_fit_status = .rc_bind_pando_field(results, "condition_fit_status"),
    pando_network_index = if (is.null(condition_result)) data.frame() else
      condition_result$pando_network_index,
    pando_fit_diagnostics = if (is.null(condition_result)) data.frame() else
      condition_result$pando_fit_diagnostics,
    pando_execution_summary = if (is.null(condition_result)) list() else
      condition_result$pando_execution_summary,
    tf_peak_gene_universal = .rc_bind_pando_field(
      results, "tf_peak_gene_universal"
    ),
    tf_peak_gene_condition_all = .rc_bind_pando_field(
      results, "tf_peak_gene_condition_all"
    ),
    tf_peak_gene_condition = .rc_bind_pando_field(
      results, "tf_peak_gene_condition"
    ),
    tf_peak_gene_condition_effect_all = if (is.null(condition_result)) {
      data.frame()
    } else condition_result$tf_peak_gene_condition_effect_all,
    tf_peak_gene_condition_effect = if (is.null(condition_result)) {
      data.frame()
    } else condition_result$tf_peak_gene_condition_effect,
    paired_cell_ids = unique(as.character(paired_meta$cell_id)),
    paired_cell_metadata = paired_meta,
    normalization_policy = list(
      rna = "global single-cell normalized RNA",
      atac = "cell-type-shared TF-IDF",
      cell_type_routing =
        "condition GRN for at least two retained conditions; standard Pando otherwise",
      condition_effect_filter = "estimable and BH adjusted P below 0.05",
      standard_edge_filter =
        "adjusted P below 0.05 and absolute estimate above 0.01",
      projection =
        "paired-cell TF-by-ATAC before exact SuperCell aggregation"
    ),
    group_cols = c(condition_col, celltype_col)
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(answer$condition_fit_status,
    file.path(outdir, "pando_group_status.tsv.gz"))
  .rc_write_tsv_gz(answer$tf_peak_gene_condition_all,
    file.path(outdir, "pando_tf_peak_gene_all.tsv.gz"))
  .rc_write_tsv_gz(answer$tf_peak_gene_condition,
    file.path(outdir, "pando_tf_peak_gene_active.tsv.gz"))
  saveRDS(answer, file.path(outdir, "single_cell_grn.rds"))
  answer
}

.rc_fit_pando_by_celltype_route <- function(
    object, gem, outdir, genome, pfm, species, condition_col, celltype_col,
    condition_types, standard_types, rna_assay, atac_assay,
    extra_args, condition_infer_args, standard_infer_args,
    parallel, BPPARAM, progress_monitor) {
  base <- list(
    object = object, gem = gem, outdir = outdir, genome = genome,
    pfm = pfm, species = species, condition_col = condition_col,
    celltype_col = celltype_col, rna_assay = rna_assay,
    atac_assay = atac_assay
  )
  build_args <- function(defaults) {
    c(defaults[setdiff(names(defaults), names(extra_args))], extra_args)
  }
  condition_result <- NULL
  if (length(condition_types)) {
    args <- build_args(base)
    args$outdir <- file.path(outdir, "condition")
    args$cell_type <- condition_types
    args$pando_infer_args <- condition_infer_args
    args$BPPARAM <- if (isTRUE(parallel)) BPPARAM else FALSE
    args$progress_monitor <- progress_monitor
    condition_result <- do.call(.rc_fit_condition_grns_by_cell_type, args)
  }
  standard_results <- list()
  for (type in standard_types) {
    cells <- rownames(object@meta.data)[
      as.character(object@meta.data[[celltype_col]]) == type
    ]
    one <- subset(object, cells = cells)
    levels <- unique(as.character(one@meta.data[[condition_col]]))
    if (length(levels) != 1L) {
      stop("Standard Pando fallback requires one retained condition in cell type `",
           type, "`.", call. = FALSE)
    }
    args <- build_args(base)
    args$object <- one
    args$outdir <- file.path(outdir, "standard", .rc_safe_path_component(type))
    args$cell_type <- type
    args$pando_infer_args <- standard_infer_args
    args$parallel <- isTRUE(parallel)
    args$progress_monitor <- progress_monitor
    standard_results[[type]] <- do.call(
      .rc_fit_standard_pando_by_cell_type, args
    )
  }
  .rc_merge_pando_results(
    condition_result = condition_result,
    standard_results = standard_results,
    condition_types = condition_types,
    standard_types = standard_types,
    condition_col = condition_col,
    celltype_col = celltype_col,
    outdir = outdir
  )
}

.rc_overlay_projection <- function(target, incoming, label) {
  if (!identical(dimnames(target), dimnames(incoming))) {
    stop(label, " projection layout is incompatible.", call. = FALSE)
  }
  occupied <- is.finite(target)
  supplied <- is.finite(incoming)
  overlap <- occupied & supplied
  if (any(overlap) &&
      any(abs(target[overlap] - incoming[overlap]) > 1e-10)) {
    stop(label, " projection overlaps another cell-type route.", call. = FALSE)
  }
  target[supplied] <- incoming[supplied]
  target
}

.rc_project_pando_by_celltype <- function(
    grn_result, membership, unit_meta, genes,
    condition_col, celltype_col, rna_assay, atac_assay) {
  template <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  projection <- reliability <- template
  coverage <- list()
  origins <- schemas <- projection_names <- policies <- character()
  if (length(grn_result$condition_grn_fits)) {
    part <- .rc_condition_pando_projection(
      grn_result, membership, unit_meta, genes
    )
    projection <- .rc_overlay_projection(
      projection, part$projection, "Condition-Pando"
    )
    reliability <- .rc_overlay_projection(
      reliability, part$reliability, "Condition-Pando reliability"
    )
    coverage[[length(coverage) + 1L]] <- part$coverage
    origins <- c(origins, part$origin)
    schemas <- c(schemas, part$pando_schema)
    projection_names <- c(projection_names, part$projection_name)
    policies <- c(policies, part$nonestimable_policy)
  }
  if (length(grn_result$standard_pando_objects)) {
    standard <- .rc_standard_pando_projection(
      grn_result, membership, unit_meta, condition_col, celltype_col,
      rna_assay = rna_assay, atac_assay = atac_assay,
      target_genes = genes
    )
    projection <- .rc_overlay_projection(
      projection, standard$projection, "Standard-Pando"
    )
    reliability <- .rc_overlay_projection(
      reliability, standard$reliability, "Standard-Pando reliability"
    )
    coverage[[length(coverage) + 1L]] <- standard$coverage
    origins <- c(origins, standard$projection_origin)
    schemas <- c(schemas, "standard_pando_network")
    projection_names <- c(projection_names, "standard_pando_full_fit")
    policies <- c(policies, "not_applicable_standard_pando")
  }
  if (!length(origins)) stop("No Pando projection route is available.", call. = FALSE)
  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = paste(unique(origins), collapse = ";"),
    pando_schema = paste(unique(schemas), collapse = ";"),
    projection_name = paste(unique(projection_names), collapse = ";"),
    nonestimable_policy = paste(unique(policies), collapse = ";"),
    condition_coefficients_calculated =
      length(grn_result$condition_grn_fits) > 0L,
    cell_type_analysis_mode = grn_result$cell_type_analysis_mode
  )
}
'''
Path("R/mixed_pando.R").write_text(mixed_file)

desc = Path("DESCRIPTION")
text = desc.read_text()
needle = "    'standard_pando.R'\n"
if needle not in text:
    raise RuntimeError("DESCRIPTION standard_pando Collate entry not found")
text = text.replace(needle, needle + "    'mixed_pando.R'\n", 1)
desc.write_text(text)

standard_path = Path("R/standard_pando.R")
text = standard_path.read_text()
old = '''      min_model_rsq = min_model_rsq,
      min_abs_estimate = max(
        .rc_standard_pando_min_abs_fixed,
        as.numeric(min_abs_estimate)
      ),
      padj_threshold = .rc_standard_pando_padj_fixed
'''
new = '''      absolute_estimate_threshold = .rc_standard_pando_min_abs_fixed,
      padj_threshold = .rc_standard_pando_padj_fixed
'''
if old not in text:
    raise RuntimeError("standard normalization-policy block not found")
standard_path.write_text(text.replace(old, new, 1))

condition_path = Path("R/condition_grn_contract.R")
text = condition_path.read_text()
text = text.replace(
'''        predictive_oof_available = FALSE,
        oof_validation_level = "not_applicable_fixed_dictionary_glm",
''', "", 1)
condition_path.write_text(text)

step_path = Path("R/step_grn_common_dictionary.R")
text = step_path.read_text()
text = text.replace(
'''  cell_type <- cell_set$retained_cell_types
''',
'''  cell_type <- cell_set$retained_cell_types
  condition_types <- cell_set$condition_pando_cell_types
  standard_types <- cell_set$standard_pando_cell_types
''', 1)
text = text.replace(
'''    skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
''',
'''    condition_pando_cell_types = condition_types,
    standard_pando_cell_types = standard_types,
''', 1)
text = text.replace(
'''  if (!identical(design$analysis_mode, cell_set$analysis_mode)) {
    stop(
      "Stage 1 cell filtering and condition design resolved different analysis modes.",
      call. = FALSE
    )
  }

''', "", 1)
branch_start = text.find('  if (identical(design$analysis_mode, "condition_grn")) {')
if branch_start < 0:
    raise RuntimeError("Stage 1 routing branch start not found")
branch_end_marker = "\n\n  grn_result$analysis_mode <- design$analysis_mode"
branch_end = text.find(branch_end_marker, branch_start)
if branch_end < 0:
    raise RuntimeError("Stage 1 routing branch end not found")
new_branch = r'''  condition_infer_args <- infer_args
  if (length(condition_types)) {
    allowed_infer_args <- c(
      "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
      "rank_action", "min_residual_df", "rna_layer", "peak_layer",
      "peak_value_type", "verbose"
    )
    unknown_infer_args <- setdiff(names(infer_args), allowed_infer_args)
    if (length(unknown_infer_args)) {
      stop("Unsupported `pando_infer_args` in condition mode: ",
           paste(unknown_infer_args, collapse = ", "), call. = FALSE)
    }
    condition_infer_args <- utils::modifyList(list(
      tf_cor = 0.1, peak_cor = 0, adjust_method = "BH",
      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L,
      verbose = .rc_progress_enabled(progress)
    ), infer_args)
    if (!identical(
          toupper(as.character(condition_infer_args$adjust_method)), "BH"
        ) ||
        !isTRUE(all.equal(
          as.numeric(condition_infer_args$padj_threshold), 0.05
        ))) {
      stop("Canonical RegCompass condition effects require BH padj < 0.05.",
           call. = FALSE)
    }
  }
  standard_infer_args <- if (length(condition_types)) {
    .rc_standard_pando_runtime_args(condition_infer_args)
  } else {
    infer_args
  }
  standard_infer_args$verbose <- standard_infer_args$verbose %||%
    .rc_progress_enabled(progress)
  grn_result <- .rc_with_step_diagnostics(
    .rc_fit_pando_by_celltype_route(
      object = object, gem = gem, outdir = outdir, genome = genome,
      pfm = pfm, species = species,
      condition_col = effective_condition_col,
      celltype_col = celltype_col,
      condition_types = condition_types,
      standard_types = standard_types,
      rna_assay = rna_assay, atac_assay = atac_assay,
      extra_args = extra_args,
      condition_infer_args = condition_infer_args,
      standard_infer_args = standard_infer_args,
      parallel = parallel, BPPARAM = BPPARAM,
      progress_monitor = monitor
    ),
    monitor
  )'''
text = text[:branch_start] + new_branch + text[branch_end:]
text = text.replace(
'''  grn_result$analysis_mode <- design$analysis_mode
''',
'''  grn_result$analysis_mode <- cell_set$analysis_mode
''', 1)
text = text.replace(
'''  grn_result$fallback_reason <- design$fallback_reason
''',
'''  grn_result$fallback_reason <- if (identical(
    cell_set$analysis_mode, "mixed_pando"
  )) "cell_type_specific_condition_count" else design$fallback_reason
''', 1)
text = text.replace(
'''      skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
''',
'''      condition_pando_cell_types = condition_types,
      standard_pando_cell_types = standard_types,
''', 1)
text = text.replace(
'''      analysis_mode = design$analysis_mode,
      fallback_reason = design$fallback_reason,
''',
'''      analysis_mode = cell_set$analysis_mode,
      fallback_reason = grn_result$fallback_reason,
      cell_type_analysis_mode = grn_result$cell_type_analysis_mode,
''', 1)
text = text.replace(
'''      pando_args = c(extra_args, list(pando_infer_args = infer_args)),
''',
'''      pando_args = c(extra_args, list(
        pando_infer_args = condition_infer_args
      )),
''', 1)
step_path.write_text(text)

layer1_path = Path("R/layer1_regulatory_support.R")
text = layer1_path.read_text()
start = text.find('  mode <- grn_result$analysis_mode %||% "condition_grn"')
end = text.find("\n  calibration <- .rc_projection_scale(", start)
if start < 0 or end < 0:
    raise RuntimeError("Layer 1 projection dispatch block not found")
replacement = r'''  mode <- grn_result$analysis_mode %||% "condition_grn"
  projection <- .rc_project_pando_by_celltype(
    grn_result = grn_result,
    membership = membership,
    unit_meta = unit_meta,
    genes = genes,
    condition_col = condition_col,
    celltype_col = celltype_col,
    rna_assay = grn_result$rna_assay %||% "RNA",
    atac_assay = grn_result$atac_assay %||% "ATAC"
  )'''
text = text[:start] + replacement + text[end:]
text = text.replace(
'''      condition_coefficients_calculated = identical(mode, "condition_grn"),
''',
'''      condition_coefficients_calculated =
        isTRUE(projection$condition_coefficients_calculated),
      cell_type_analysis_mode = projection$cell_type_analysis_mode,
''', 1)
layer1_path.write_text(text)

stepwise = Path("R/stepwise_workflow.R")
text = stepwise.read_text()
text = text.replace(
'''      skipped_condition_cell_types = cell_set$skipped_condition_cell_types,
''',
'''      condition_pando_cell_types = cell_set$condition_pando_cell_types,
      standard_pando_cell_types = cell_set$standard_pando_cell_types,
''', 1)
text = text.replace(
'''  effective_condition_col <- design$condition_col
''',
'''  effective_condition_col <- design$condition_col
  resolved_analysis_mode <- if (is.null(grn)) {
    contract$analysis_mode
  } else {
    grn$params$analysis_mode
  }
  resolved_fallback_reason <- if (is.null(grn)) {
    design$fallback_reason
  } else {
    grn$params$fallback_reason
  }
''', 1)
text = text.replace(
'''      analysis_mode = design$analysis_mode,
      fallback_reason = design$fallback_reason,
''',
'''      analysis_mode = resolved_analysis_mode,
      fallback_reason = resolved_fallback_reason,
      cell_type_analysis_mode = if (is.null(grn)) NULL else
        grn$params$cell_type_analysis_mode,
''', 1)
stepwise.write_text(text)

results = Path("R/step_results.R")
text = results.read_text()
text = text.replace(
'''    condition_coefficients_calculated = identical(mode, "condition_grn"),
''',
'''    condition_coefficients_calculated =
      isTRUE(grn$grn_result$condition_coefficients_calculated),
''', 1)
old_design = '''      pando_design = if (identical(mode, "condition_grn")) {
        paste(
          "pooled and condition-specific candidate discovery, exact",
          "TF-peak-target union, and one unscaled fixed-dictionary Gaussian",
          "identity GLM per condition"
        )
      } else {
        "original Pando infer_grn per broad cell type; no condition coefficients"
      },
'''
new_design = '''      pando_design =
        "cell-type-specific routing: common-dictionary condition GRN for at least two retained conditions, standard Pando otherwise",
      cell_type_analysis_mode = grn$grn_result$cell_type_analysis_mode,
'''
if old_design not in text:
    raise RuntimeError("result pando_design block not found")
results.write_text(text.replace(old_design, new_design, 1))

reg = Path("R/regcompass.R")
text = reg.read_text()
text = text.replace(
'''  result$params$condition_coefficients_calculated <-
    identical(step1$params$analysis_mode, "condition_grn")
''',
'''  result$params$condition_coefficients_calculated <-
    isTRUE(step1$grn_result$condition_coefficients_calculated)
  result$params$cell_type_analysis_mode <-
    step1$grn_result$cell_type_analysis_mode
''', 1)
reg.write_text(text)

replace_once(
    "R/step_layer1.R",
    "#' Condition mode uses outer-heldout condition projections. Standard mode uses\n",
    "#' Each cell type uses its Stage 1 Pando route. Standard mode uses\n",
    "Layer 1 roxygen"
)

source = "\n".join(p.read_text() for p in Path("R").glob("*.R"))
for token in (
    "condition_full_oof", "common_oof", "condition_unique_oof",
    'network_name = "regcompass_condition_grn"',
    "min_model_rsq", "min_abs_estimate",
    "predictive_oof_available", "oof_validation_level"
):
    if token in source:
        raise RuntimeError(f"retired token remains in R source: {token}")
