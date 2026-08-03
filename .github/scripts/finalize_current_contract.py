from __future__ import annotations

from pathlib import Path
import re

ROOT = Path('.')
PANDO_SHA = '05246162ab5639fb12407b2b329de5149d9660a4'


def function_span(text: str, name: str) -> tuple[int, int]:
    pattern = re.compile(r'(?m)^' + re.escape(name) + r'\s*<-\s*function\b')
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(f'{name}: expected one definition, found {len(matches)}')
    start = matches[0].start()
    opening = text.find('{', matches[0].end())
    if opening < 0:
        raise RuntimeError(f'{name}: opening brace not found')
    depth = 0
    quote = None
    escaped = False
    comment = False
    for index in range(opening, len(text)):
        char = text[index]
        if comment:
            if char == '\n':
                comment = False
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == '\\':
                escaped = True
            elif char == quote:
                quote = None
            continue
        if char == '#':
            comment = True
        elif char in ('"', "'"):
            quote = char
        elif char == '{':
            depth += 1
        elif char == '}':
            depth -= 1
            if depth == 0:
                return start, index + 1
    raise RuntimeError(f'{name}: closing brace not found')


def replace_function(path: str, name: str, replacement: str) -> None:
    file = Path(path)
    text = file.read_text()
    start, end = function_span(text, name)
    file.write_text(text[:start] + replacement.rstrip() + text[end:])


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f'{label}: expected one match, found {count}')
    return text.replace(old, new, 1)


# ---------------------------------------------------------------------------
# Layer 1: one projection, one modifier, one multiome reaction matrix.
# ---------------------------------------------------------------------------
replace_function(
    'R/step_layer1_common_dictionary.R',
    '.rc_condition_pando_projection',
    r'''.rc_condition_pando_projection <- function(
    grn_result, membership, unit_meta, genes) {
  projection <- matrix(
    NA_real_, length(genes), nrow(unit_meta),
    dimnames = list(genes, unit_meta$unit_id)
  )
  reliability <- projection
  coverage <- list()

  for (fit in grn_result$condition_grn_fits) {
    .rc_require_pando_condition_grn_fit(fit)
    cell_projection <- Pando::project_condition_grn_cells(
      object = grn_result$pando_grn_data,
      fit = fit,
      targets = genes,
      significant_only = TRUE,
      return_edge_contributions = FALSE
    )
    aggregated <- Pando::aggregate_condition_grn_projection(
      cell_projection, membership, group_col = "metacell_id"
    )
    score <- as.matrix(aggregated$gene_score)
    rownames(score) <- tolower(rownames(score))
    targets <- intersect(rownames(score), rownames(projection))
    units <- intersect(colnames(score), colnames(projection))
    projection[targets, units] <- score[targets, units, drop = FALSE]

    coefficient <- as.data.frame(fit$coefficients, stringsAsFactors = FALSE)
    condition_col <- as.character(fit$condition_col)[[1L]]
    celltype_col <- as.character(fit$cell_type_col)[[1L]]
    if (!all(c(condition_col, celltype_col) %in% colnames(unit_meta))) {
      stop("Metacell metadata lack fitted condition or cell-type columns.",
           call. = FALSE)
    }
    for (condition in fit$condition_levels) {
      selected_units <- unit_meta$unit_id[
        as.character(unit_meta[[condition_col]]) == condition &
          as.character(unit_meta[[celltype_col]]) == fit$cell_type
      ]
      significant_edges <- coefficient[
        as.character(coefficient$condition) == condition &
          coefficient$estimable %in% TRUE &
          is.finite(as.numeric(coefficient$padj)) &
          as.numeric(coefficient$padj) < 0.05,
        , drop = FALSE
      ]
      reliable_targets <- intersect(
        tolower(unique(as.character(significant_edges$target))),
        rownames(reliability)
      )
      if (length(reliable_targets) && length(selected_units)) {
        reliability[reliable_targets, selected_units] <- 1
      }
    }
    coverage[[length(coverage) + 1L]] <- data.frame(
      cell_type = fit$cell_type,
      condition = fit$condition_levels,
      n_dictionary_edges = nrow(fit$edge_dictionary),
      n_significant_edges = vapply(fit$condition_levels, function(condition) {
        sum(
          coefficient$condition == condition &
            coefficient$estimable %in% TRUE &
            is.finite(as.numeric(coefficient$padj)) &
            as.numeric(coefficient$padj) < 0.05,
          na.rm = TRUE
        )
      }, integer(1)),
      padj_threshold = 0.05,
      projection_effect = "penalty_effect",
      stringsAsFactors = FALSE
    )
  }

  list(
    projection = projection,
    reliability = reliability,
    coverage = .rc_bind_frames_fill(coverage),
    origin = "paired_cell_fixed_dictionary_glm_padj_filtered",
    pando_schema = .RC_PANDO_CONDITION_GRN_FIT_SCHEMA,
    projection_name = "padj_filtered_fixed_dictionary_condition_glm",
    nonestimable_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution"
  )
}'''
)

replace_function(
    'R/layer1_regulatory_support.R',
    '.rc_latent_metacell_expression',
    r'''.rc_latent_metacell_expression <- function(
    counts, library_size, mu_min = 0.1, cell_type) {
  counts <- as.matrix(counts)
  library_size <- as.numeric(library_size)
  if (!is.numeric(counts) || is.null(rownames(counts)) ||
      is.null(colnames(counts)) || length(library_size) != ncol(counts) ||
      any(!is.finite(library_size)) || any(library_size <= 0) ||
      any(!is.finite(counts)) || any(counts < 0)) {
    stop("Latent RNA estimation requires aligned non-negative counts and depth.",
         call. = FALSE)
  }
  if (any(abs(counts - round(counts)) >
      sqrt(.Machine$double.eps) * pmax(1, abs(counts)))) {
    stop("Latent RNA estimation requires raw integer-like counts.",
         call. = FALSE)
  }
  if (!is.null(names(cell_type))) {
    cell_type <- as.character(cell_type[colnames(counts)])
  } else {
    cell_type <- as.character(cell_type)
  }
  if (length(cell_type) != ncol(counts) || anyNA(cell_type)) {
    stop("Cell-type labels do not align to metacell RNA counts.",
         call. = FALSE)
  }
  exposure <- library_size / 1e6
  observed_cpm <- sweep(counts, 2L, exposure, "/")
  types <- unique(cell_type)
  prior_shape <- prior_rate <- prior_mean <- prior_variance <- matrix(
    NA_real_, nrow(counts), length(types),
    dimnames = list(rownames(counts), types)
  )
  for (value in types) {
    x <- observed_cpm[, cell_type == value, drop = FALSE]
    mean_value <- rowMeans(x)
    variance_value <- if (ncol(x) > 1L) {
      apply(x, 1L, stats::var)
    } else {
      rep(NA_real_, nrow(x))
    }
    invalid <- !is.finite(variance_value) | variance_value <= 0
    variance_value[invalid] <- pmax(mean_value[invalid]^2 / 1e6, 1e-12)
    shape_value <- pmin(
      1e6, pmax(0.1, mean_value^2 / pmax(variance_value, 1e-12))
    )
    rate_value <- shape_value / pmax(mean_value, 1e-8)
    prior_mean[, value] <- mean_value
    prior_variance[, value] <- variance_value
    prior_shape[, value] <- shape_value
    prior_rate[, value] <- rate_value
  }
  shape_unit <- prior_shape[, cell_type, drop = FALSE]
  rate_unit <- prior_rate[, cell_type, drop = FALSE]
  colnames(shape_unit) <- colnames(rate_unit) <- colnames(counts)
  posterior_shape <- observed_cpm + shape_unit
  posterior_rate <- rate_unit + 1
  latent_cpm <- posterior_shape / posterior_rate
  positive <- matrix(
    stats::pgamma(
      rep(mu_min, length(posterior_shape)),
      shape = as.numeric(posterior_shape),
      rate = as.numeric(posterior_rate),
      lower.tail = FALSE
    ),
    nrow = nrow(counts), dimnames = dimnames(counts)
  )
  observed_zero <- counts == 0
  zero_class <- matrix(
    "observed_positive", nrow(counts), ncol(counts),
    dimnames = dimnames(counts)
  )
  zero_class[observed_zero] <- "observed_zero_continuous_eb"
  observation_weight <- 1 / posterior_rate
  prior_weight <- rate_unit / posterior_rate
  list(
    latent_cpm = latent_cpm,
    latent_log_expression = log1p(latent_cpm),
    posterior_positive_probability = positive,
    posterior_zero_probability = 1 - positive,
    zero_class = zero_class,
    observed_zero = observed_zero,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_shape = prior_shape,
    prior_rate = prior_rate,
    prior_weight = prior_weight,
    observation_weight = observation_weight,
    model = "normalized_unit_gamma_empirical_bayes_by_cell_type_v4",
    prior_estimation_scope = "gene_by_cell_type",
    posterior_update_scope =
      "one_normalized_metacell_unit_independent_of_library_size"
  )
}'''
)

replace_function(
    'R/layer1_regulatory_support.R',
    '.rc_cell_first_projection_layer1',
    r'''.rc_cell_first_projection_layer1 <- function(
    grn_result, metacell_object, membership, metacell_meta, gem,
    condition_col, celltype_col, rna_assay,
    gpr_and_method = "min", gene_half_saturation = 1,
    parallel = TRUE, BPPARAM = NULL) {
  if (!is.data.frame(membership) ||
      !all(c("cell_id", "metacell_id") %in% colnames(membership)) ||
      anyDuplicated(membership$cell_id)) {
    stop("SuperCell membership must map every cell exactly once.",
         call. = FALSE)
  }
  id_col <- if ("metacell_id" %in% colnames(metacell_meta)) {
    "metacell_id"
  } else {
    "pool_id"
  }
  unit_meta <- as.data.frame(metacell_meta)
  unit_meta$unit_id <- as.character(unit_meta[[id_col]])
  unit_meta$pool_id <- unit_meta$unit_id
  units <- colnames(metacell_object)
  unit_meta <- unit_meta[match(units, unit_meta$unit_id), , drop = FALSE]
  if (anyNA(unit_meta$unit_id)) {
    stop("Metacell metadata do not align to the metacell object.",
         call. = FALSE)
  }
  parsed <- rc_parse_gpr_table(gem$gpr_table)
  gpr_genes <- unique(tolower(unlist(parsed, use.names = FALSE)))
  counts <- .rc_get_assay_counts(metacell_object, rna_assay)
  library_size <- Matrix::colSums(counts)
  rna_counts <- counts[
    tolower(rownames(counts)) %in% gpr_genes, units, drop = FALSE
  ]
  rownames(rna_counts) <- tolower(rownames(rna_counts))
  if (anyDuplicated(rownames(rna_counts))) {
    stop("Duplicated GPR genes after case normalization.", call. = FALSE)
  }
  cell_type <- stats::setNames(
    as.character(unit_meta[[celltype_col]]), unit_meta$unit_id
  )
  latent <- .rc_latent_metacell_expression(
    rna_counts, library_size[units], cell_type = cell_type
  )
  genes <- rownames(latent$latent_log_expression)
  mode <- grn_result$analysis_mode %||% "condition_grn"
  projection <- if (identical(mode, "standard_pando")) {
    standard <- .rc_standard_pando_projection(
      grn_result, membership, unit_meta, condition_col, celltype_col,
      rna_assay = grn_result$rna_assay %||% "RNA",
      atac_assay = grn_result$atac_assay %||% "ATAC",
      target_genes = genes
    )
    list(
      projection = standard$projection,
      reliability = standard$reliability,
      coverage = standard$coverage,
      origin = standard$projection_origin,
      pando_schema = "standard_pando_network",
      projection_name = "standard_pando_full_fit",
      nonestimable_policy = "not_applicable_standard_pando"
    )
  } else {
    .rc_condition_pando_projection(
      grn_result, membership, unit_meta, genes
    )
  }
  calibration <- .rc_projection_scale(
    projection$projection, unit_meta, celltype_col
  )
  modifier <- .rc_scaled_regulatory_modifier(
    projection$projection, projection$reliability, calibration$scale
  )
  gene_support_rna <- rc_gene_score(
    latent$latent_log_expression,
    mode = "absolute",
    half_saturation = gene_half_saturation
  )
  gene_support_multiome <- .rc_integrate_regulatory_support(
    gene_support_rna, modifier, alpha = 1
  )
  reaction_rna <- rc_reaction_capacity(
    parsed, gene_support_rna, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  reaction_multiome <- rc_reaction_capacity(
    parsed, gene_support_multiome, promiscuity_mode = "none",
    and_method = gpr_and_method, or_method = "sum",
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  support_fraction <- .rc_gpr_best_group_fraction(
    parsed, is.finite(modifier)
  )
  fallback <- !is.finite(modifier)
  list(
    schema_version = "regcompass_regulatory_layer1_v4",
    analysis_mode = mode,
    reaction_expression = reaction_multiome,
    reaction_expression_rna_only = reaction_rna,
    reaction_expression_available = is.finite(reaction_multiome),
    rna_metacell_latent_log_expression = latent$latent_log_expression,
    rna_metacell_latent_cpm = latent$latent_cpm,
    posterior_positive_probability = latent$posterior_positive_probability,
    posterior_zero_probability = latent$posterior_zero_probability,
    rna_zero_class = latent$zero_class,
    eb_prior_weight = latent$prior_weight,
    eb_observation_weight = latent$observation_weight,
    gene_support_rna = gene_support_rna,
    gene_support_multiome = gene_support_multiome,
    gene_projection = projection$projection,
    gene_projection_scale = calibration$scale,
    gene_regulatory_reliability = projection$reliability,
    gene_regulatory_reliability_available =
      is.finite(projection$reliability),
    gene_regulatory_available = is.finite(projection$projection),
    gene_regulatory_modifier = modifier,
    projection_coverage = projection$coverage,
    projection_calibration = calibration$diagnostics,
    reaction_regulatory_support_fraction = support_fraction,
    parsed_gpr = parsed,
    gpr_diagnostics = rc_gpr_diagnostics(parsed, genes),
    unit_meta = unit_meta,
    metacell_meta = unit_meta,
    layer1_unit = "native_SuperCell_metacell",
    regulatory_fallback = list(
      policy = "rna_only_for_nonfinite_pando_modifier",
      neutral_modifier = 0,
      gene_metacell_mask = fallback,
      n_fallback = sum(fallback),
      fallback_fraction = mean(fallback)
    ),
    depth_diagnostics = list(
      rna_library_size = stats::setNames(
        as.numeric(library_size[units]), units
      ),
      latent_expression_model = latent$model,
      prior_estimation_scope = latent$prior_estimation_scope,
      posterior_update_scope = latent$posterior_update_scope,
      eb_weight_library_size_dependent = FALSE
    ),
    zero_diagnostics = list(
      observed_zero_fraction = rowMeans(latent$observed_zero),
      mean_posterior_zero_probability = rowMeans(
        latent$posterior_zero_probability
      ),
      maximum_posterior_zero_probability = apply(
        latent$posterior_zero_probability, 1L, max
      )
    ),
    capacity_params = list(
      regulatory_odds_budget = 2,
      gene_half_saturation = gene_half_saturation,
      regulatory_mode = mode,
      link_function = "tanh(G/shared_scale)",
      promiscuity_mode = "none",
      and_method = gpr_and_method,
      or_method = "sum",
      parallel = parallel
    ),
    projection_provenance = list(
      analysis_mode = mode,
      pando_schema = projection$pando_schema,
      projection_origin = projection$origin,
      projection_used_for_penalty = TRUE,
      projection_name = projection$projection_name,
      condition_coefficients_calculated = identical(mode, "condition_grn"),
      supercell_membership = "membership_table(cell_id, metacell_id)",
      unavailable_target_policy = "rna_only_neutral_modifier_fallback",
      nonestimable_edge_policy = projection$nonestimable_policy
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    statistical_unit = "metacell",
    metacell_statistical_inference = TRUE,
    biological_replicate_inference = FALSE
  )
}'''
)

replace_function(
    'R/step_layer1.R',
    'rc_regcompass_step_layer1',
    r'''rc_regcompass_step_layer1 <- function(
    grn, metacells, meta_modules, gem, outdir,
    gpr_and_method = c("min", "median", "mean"),
    gene_half_saturation = getOption("RegCompassR.cpm_half_saturation", 1),
    parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer1", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gpr_and_method <- match.arg(gpr_and_method)
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!identical(metacells$params, meta_modules$workflow_params)) {
    stop("Metacell and meta-module stages use different workflow settings.",
         call. = FALSE)
  }
  if (!identical(grn$params$analysis_mode, metacells$params$analysis_mode)) {
    stop("GRN and metacell stages resolved different analysis modes.",
         call. = FALSE)
  }
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_require_stage_gem(grn, gem, "grn")
  params <- metacells$params
  layer1 <- .rc_cell_first_projection_layer1(
    grn_result = grn$grn_result,
    metacell_object = metacells$metacell_object,
    membership = metacells$pooled$membership,
    gem = gem,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    rna_assay = params$rna_assay,
    gpr_and_method = gpr_and_method,
    gene_half_saturation = gene_half_saturation,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  layer1$workflow_params <- params
  layer1$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(layer1) <- c("regcompass_layer1_step", "list")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer1 <- .rc_step_monitor_finish(layer1, monitor)
  saveRDS(layer1, file.path(outdir, "step_layer1.rds"))
  layer1
}'''
)

# ---------------------------------------------------------------------------
# Layer 2: current penalty plus one RNA-only interpretation control.
# ---------------------------------------------------------------------------
replace_function(
    'R/step_layer2.R',
    '.rc_layer2_comparison_table',
    r'''.rc_layer2_comparison_table <- function(
    layer2, layer1, condition_col, celltype_col) {
  penalty <- layer2$penalty
  row_meta <- rc_parse_microcompass_row_id(rownames(penalty))
  unit_meta <- layer2$unit_meta
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else {
    as.character(unit_meta$pool_id)
  }
  unit_meta <- unit_meta[match(colnames(penalty), unit_id), , drop = FALSE]
  unit_celltype <- as.character(unit_meta[[celltype_col]])
  grid <- expand.grid(
    row_index = seq_len(nrow(penalty)),
    unit_index = seq_len(ncol(penalty)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scoped_celltype <- row_meta$cell_type[grid$row_index]
  matching_scope <- is.na(scoped_celltype) |
    scoped_celltype == unit_celltype[grid$unit_index]
  grid <- grid[matching_scope, , drop = FALSE]
  if (!nrow(grid)) {
    stop("No Layer 2 rows match their cell-type units.", call. = FALSE)
  }
  reaction <- row_meta$reaction_id[grid$row_index]
  unit <- colnames(penalty)[grid$unit_index]
  index_matrix <- function(x) {
    as.numeric(x[cbind(grid$row_index, grid$unit_index)])
  }
  reaction_unit_matrix <- function(x) {
    if (!is.numeric(x) || is.null(dim(x)) ||
        !all(unique(reaction) %in% rownames(x)) ||
        !identical(colnames(x), colnames(penalty))) {
      stop("Layer 1 reaction diagnostics are not aligned.", call. = FALSE)
    }
    as.numeric(x[cbind(match(reaction, rownames(x)), grid$unit_index)])
  }
  omega <- layer2$params$omega
  vmax <- index_matrix(layer2$vmax)
  primary <- index_matrix(layer2$penalty)
  rna_only <- index_matrix(layer2$penalty_rna_only)
  normalized <- primary / (omega * vmax)
  normalized[!is.finite(primary) | !is.finite(vmax) | vmax <= 0] <- NA_real_
  data.frame(
    row_id = rownames(penalty)[grid$row_index],
    reaction_id = reaction,
    direction = row_meta$target_direction[grid$row_index],
    medium = row_meta$medium_scenario[grid$row_index],
    cell_type = unit_celltype[grid$unit_index],
    condition = as.character(unit_meta[[condition_col]][grid$unit_index]),
    metacell_id = unit,
    penalty = primary,
    penalty_rna_only = rna_only,
    regulatory_penalty_delta = primary - rna_only,
    penalty_per_target_flux = normalized,
    vmax = vmax,
    penalty_available = is.finite(primary),
    regulatory_support_fraction = reaction_unit_matrix(
      layer1$reaction_regulatory_support_fraction
    ),
    inference_class = "metacell_statistical_unit_within_dataset",
    comparability_class =
      "same_celltype_conditions_on_one_celltype_medium_union_gem",
    stringsAsFactors = FALSE
  )
}'''
)

layer2_path = Path('R/step_layer2.R')
layer2_text = layer2_path.read_text()
layer2_text = layer2_text.replace(
    "# Layer 2 shared-structure validation and condition-full scoring.",
    "# Layer 2 shared-structure validation and current penalty scoring.",
    1
)
layer2_text = layer2_text.replace(
    "#' share a model only when their cell type matches. Historical `_oof`, `common`,\n#' and `condition_unique` fields remain compatibility aliases; RNA-only scoring\n#' is an interpretation control. The `full_gem` route uses a separate engine.\n",
    "#' share a model only when their cell type matches. RNA-only scoring is an\n#' interpretation control that reuses the exact structural model cache. The\n#' `full_gem` route uses a separate engine.\n",
    1
)
old_tail = '''  common <- run_control(
    "reaction_expression_common_oof",
    "common-support compatibility route"
  )
  rna_only <- run_control(
    "reaction_expression_rna_only",
    "RNA-only control"
  )
  answer$schema_version <- "regcompass_regulatory_layer2_v2"
  answer$penalty_condition_full_oof <- answer$penalty
  answer$penalty_common_oof <- common$penalty
  answer$penalty_condition_unique_increment <-
    answer$penalty_condition_full_oof - answer$penalty_common_oof
  answer$penalty_rna_only <- rna_only$penalty
  answer$score_condition_full_oof_display_only <- answer$score
  answer$score_common_oof_display_only <- common$score
  answer$score_rna_only_display_only <- rna_only$score
  answer$comparison_contract <- list(
    primary = "penalty_condition_full_oof",
    common_component = "penalty_common_oof",
    condition_unique_increment =
      "penalty_condition_full_oof - penalty_common_oof",
    rna_control = "penalty_rna_only",
    nonestimable_edge_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    removed_guardrails = c(
      "depth_matching",
      "common_depth_restriction",
      "alpha_sensitivity",
      "zero_support_sensitivity",
      "link_saturation_propagation"
    ),
    exact_shared_structure = TRUE,
    structural_model_contract = answer$structural_model_contract,
    effect_size_basis = "penalty / (omega * vmax)",
    ecdf_effect_size_eligible = FALSE
  )
  common$shared_model_cache <- NULL
  rna_only$shared_model_cache <- NULL
  answer$comparison_paths <- list(
    common_support = common,
    rna_only = rna_only
  )
'''
new_tail = '''  rna_only <- run_control(
    "reaction_expression_rna_only",
    "RNA-only control"
  )
  answer$schema_version <- "regcompass_regulatory_layer2_v3"
  answer$penalty_rna_only <- rna_only$penalty
  answer$score_rna_only <- rna_only$score
  answer$comparison_contract <- list(
    primary = "penalty",
    rna_control = "penalty_rna_only",
    nonestimable_edge_policy =
      "coefficient_NA_and_zero_realized_penalty_contribution",
    exact_shared_structure = TRUE,
    structural_model_contract = answer$structural_model_contract,
    effect_size_basis = "penalty / (omega * vmax)",
    ecdf_effect_size_eligible = FALSE
  )
  rna_only$shared_model_cache <- NULL
  answer$comparison_paths <- list(rna_only = rna_only)
'''
layer2_text = replace_once(layer2_text, old_tail, new_tail, 'Layer 2 route cleanup')
layer2_path.write_text(layer2_text)

# ---------------------------------------------------------------------------
# Stage validators follow only the current schemas.
# ---------------------------------------------------------------------------
replace_function(
    'R/stage_contracts.R',
    '.rc_validate_layer1_stage',
    r'''.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL, argument = "layer1") {
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument, "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  mode <- layer1$analysis_mode
  if (!mode %in% c("condition_grn", "standard_pando")) {
    stop("Layer 1 analysis mode is invalid.", call. = FALSE)
  }
  required_gene_matrices <- c(
    "gene_projection", "gene_projection_scale",
    "gene_regulatory_reliability",
    "gene_regulatory_reliability_available",
    "gene_regulatory_modifier", "gene_support_rna",
    "gene_support_multiome"
  )
  reference <- layer1$gene_projection
  if (!is.numeric(reference) || is.null(dim(reference)) ||
      !identical(colnames(reference), ids)) {
    stop("Layer 1 regulatory projection is invalid.", call. = FALSE)
  }
  for (name in required_gene_matrices) {
    value <- layer1[[name]]
    if (is.null(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(reference))) {
      stop("Layer 1 `", name, "` is missing or misaligned.",
           call. = FALSE)
    }
  }
  if (!is.logical(layer1$gene_regulatory_reliability_available) ||
      anyNA(layer1$gene_regulatory_reliability_available)) {
    stop("Layer 1 reliability availability must be logical.",
         call. = FALSE)
  }
  expected_modifier <- layer1$gene_regulatory_reliability *
    tanh(reference / layer1$gene_projection_scale)
  expected_modifier[
    !is.finite(reference) |
      !is.finite(layer1$gene_regulatory_reliability)
  ] <- NA_real_
  observed <- layer1$gene_regulatory_modifier
  finite <- is.finite(expected_modifier) & is.finite(observed)
  if (any(is.finite(expected_modifier) != is.finite(observed)) ||
      any(abs(expected_modifier[finite] - observed[finite]) > 1e-10)) {
    stop("Layer 1 modifier is not reliability*tanh(projection/scale).",
         call. = FALSE)
  }
  for (name in c("reaction_expression_rna_only")) {
    value <- layer1[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(layer1$reaction_expression))) {
      stop("Layer 1 reaction matrices are misaligned.", call. = FALSE)
    }
  }
  if (!is.numeric(layer1$reaction_regulatory_support_fraction) ||
      !identical(
        dimnames(layer1$reaction_regulatory_support_fraction),
        dimnames(layer1$reaction_expression)
      )) {
    stop("Layer 1 reaction regulatory support is misaligned.",
         call. = FALSE)
  }
  retired <- c(
    "reaction_expression_condition_full_oof",
    "reaction_expression_common_oof",
    "gene_projection_condition_full_oof",
    "gene_projection_common_oof",
    "gene_projection_condition_unique_oof",
    "gene_support_condition_full_oof",
    "gene_support_common_oof",
    "gene_regulatory_modifier_condition_full_oof",
    "gene_regulatory_modifier_common_oof",
    "reaction_condition_full_support_fraction",
    "reaction_common_support_fraction"
  )
  if (any(retired %in% names(layer1))) {
    stop("Layer 1 contains retired projection routes.", call. = FALSE)
  }
  provenance <- layer1$projection_provenance
  expected <- if (identical(mode, "condition_grn")) {
    list(
      origin = "paired_cell_fixed_dictionary_glm_padj_filtered",
      projection = "padj_filtered_fixed_dictionary_condition_glm",
      nonestimable =
        "coefficient_NA_and_zero_realized_penalty_contribution",
      coefficients = TRUE
    )
  } else {
    list(
      origin = "standard_pando_full_fit",
      projection = "standard_pando_full_fit",
      nonestimable = "not_applicable_standard_pando",
      coefficients = FALSE
    )
  }
  if (!identical(layer1$schema_version, "regcompass_regulatory_layer1_v4") ||
      !is.list(provenance) ||
      !identical(provenance$analysis_mode, mode) ||
      !identical(provenance$projection_origin, expected$origin) ||
      !identical(provenance$projection_name, expected$projection) ||
      !identical(provenance$nonestimable_edge_policy, expected$nonestimable) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !identical(
        provenance$condition_coefficients_calculated,
        expected$coefficients
      ) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      )) {
    stop("Layer 1 regulatory projection provenance is incompatible.",
         call. = FALSE)
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer1, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer1, gem, argument)
  invisible(ids)
}'''
)

replace_function(
    'R/stage_contracts.R',
    '.rc_validate_layer2_stage',
    r'''.rc_validate_layer2_stage <- function(
    layer2, layer1 = NULL, workflow_params = NULL, gem = NULL,
    argument = "layer2", required_mode = NULL) {
  .rc_require_stage_class(
    layer2, "regcompass_layer2_step", argument, "rc_regcompass_step_layer2"
  )
  ids <- .rc_layer2_unit_ids(layer2)
  valid_modes <- c("meta_module_gem", "full_gem")
  if (!is.character(layer2$model_mode) || length(layer2$model_mode) != 1L ||
      is.na(layer2$model_mode) || !layer2$model_mode %in% valid_modes) {
    stop("Layer 2 model mode is invalid.", call. = FALSE)
  }
  if (!is.null(required_mode)) {
    if (!is.character(required_mode) || length(required_mode) != 1L ||
        is.na(required_mode) || !required_mode %in% valid_modes) {
      stop(
        "`required_mode` must be one of: ",
        paste(valid_modes, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!identical(layer2$model_mode, required_mode)) {
      stop(
        "Layer 2 model mode is `", layer2$model_mode,
        "`; `", required_mode, "` is required.", call. = FALSE
      )
    }
  }
  reference <- layer2$penalty
  for (name in c("penalty_rna_only", "score_rna_only")) {
    value <- layer2[[name]]
    if (!is.numeric(value) || is.null(dim(value)) ||
        !identical(dimnames(value), dimnames(reference))) {
      stop("Layer 2 route `", name, "` is missing or misaligned.",
           call. = FALSE)
    }
  }
  retired <- c(
    "penalty_condition_full_oof", "penalty_common_oof",
    "penalty_condition_unique_increment",
    "score_condition_full_oof_display_only",
    "score_common_oof_display_only",
    "score_rna_only_display_only"
  )
  if (any(retired %in% names(layer2))) {
    stop("Layer 2 contains retired penalty routes.", call. = FALSE)
  }
  contract <- layer2$comparison_contract
  required_contract <- c(
    "primary", "rna_control", "nonestimable_edge_policy",
    "exact_shared_structure", "structural_model_contract",
    "effect_size_basis", "ecdf_effect_size_eligible"
  )
  if (!identical(layer2$schema_version, "regcompass_regulatory_layer2_v3") ||
      !is.list(contract) ||
      !all(required_contract %in% names(contract)) ||
      !identical(contract$primary, "penalty") ||
      !identical(contract$rna_control, "penalty_rna_only") ||
      !isTRUE(contract$exact_shared_structure) ||
      !identical(
        contract$nonestimable_edge_policy,
        "coefficient_NA_and_zero_realized_penalty_contribution"
      )) {
    stop("Layer 2 comparison contract is incomplete.", call. = FALSE)
  }
  if (!is.null(layer1)) {
    .rc_validate_layer1_stage(layer1, argument = "layer1")
    if (!identical(ids, .rc_layer1_unit_ids(layer1))) {
      stop("Layer 1 and Layer 2 unit IDs differ.", call. = FALSE)
    }
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer2, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer2, gem, argument)
  invisible(ids)
}'''
)

# ---------------------------------------------------------------------------
# Pando contracts: remove compatibility arguments and old threshold gates.
# ---------------------------------------------------------------------------
replace_function(
    'R/condition_grn_contract.R',
    '.rc_pando_execution_summary',
    r'''.rc_pando_execution_summary <- function(diagnostics = NULL) {
  list(
    fit_engine = "two_stage_exact_edge_union_fixed_dictionary_glm",
    targets_total = if (is.data.frame(diagnostics)) {
      length(unique(as.character(diagnostics$target)))
    } else 0L,
    targets_failed = if (is.data.frame(diagnostics) &&
      "fit_status" %in% colnames(diagnostics)) {
      sum(!diagnostics$fit_status %in% c("ok", "rank_deficient"), na.rm = TRUE)
    } else 0L
  )
}'''
)

condition_path = Path('R/condition_grn_contract.R')
condition_text = condition_path.read_text()
condition_text = replace_once(
    condition_text,
    '''  fits <- Pando::condition_grn_fit(
    grn_object, network_name = "regcompass_condition_grn"
  )
''',
    '''  fits <- Pando::condition_grn_fit(grn_object)
''',
    'condition fit compatibility argument'
)
condition_text = replace_once(
    condition_text,
    '''    "padj_threshold", "adjust_method", "scale", "interaction",
    "projection_effect_column", "projection_policy"
''',
    '''    "padj_threshold", "adjust_method", "scale", "interaction",
    "projection_effect_column", "projection_policy", "rna_layer",
    "peak_layer", "peak_value_type", "preprocessing_fingerprint",
    "dictionary_preprocessing_provenance_verified"
''',
    'Pando fit provenance requirements'
)
condition_text = replace_once(
    condition_text,
    '''      !identical(toupper(as.character(fit$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(fit$padj_threshold), 0.05))) {
''',
    '''      !identical(toupper(as.character(fit$adjust_method)), "BH") ||
      !isTRUE(all.equal(as.numeric(fit$padj_threshold), 0.05)) ||
      !isTRUE(fit$dictionary_preprocessing_provenance_verified) ||
      any(!nzchar(c(
        as.character(fit$rna_layer), as.character(fit$peak_layer),
        as.character(fit$peak_value_type),
        as.character(fit$preprocessing_fingerprint)
      )))) {
''',
    'Pando fit provenance validation'
)
old_retired = '''  retired <- intersect(names(pando_infer_args), c(
    "condition_mix", "condition_weight", "alpha", "nlambda", "lambda",
    "lambda_min_ratio", "outer_nfolds", "inner_nfolds",
    "lambda_selection", "scale", "active_tol", "max_iter",
    "tol_objective", "tol_coef", "seed", "comparison_conditions"
  ))
  if (length(retired)) {
    stop(
      "Retired nested-CV condition-GRN argument(s): ",
      paste(retired, collapse = ", "),
      ". Use tf_cor, peak_cor, adjust_method, padj_threshold, rank_action and min_residual_df.",
      call. = FALSE
    )
  }
  pando_infer_args$candidate_screen <- NULL
  pando_infer_args$engine_control <- NULL
  pando_infer_args$parallel <- NULL
  pando_infer_args$verbose <- NULL
'''
new_allowed = '''  allowed_infer_args <- c(
    "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
    "rank_action", "min_residual_df", "rna_layer", "peak_layer",
    "peak_value_type"
  )
  unknown_infer_args <- setdiff(names(pando_infer_args), allowed_infer_args)
  if (length(unknown_infer_args)) {
    stop(
      "Unsupported `pando_infer_args`: ",
      paste(unknown_infer_args, collapse = ", "), call. = FALSE
    )
  }
'''
condition_text = replace_once(
    condition_text, old_retired, new_allowed,
    'condition inference argument cleanup'
)
condition_path.write_text(condition_text)

replace_function(
    'R/pando_grn.R',
    '.rc_validate_pando_evidence_filters',
    r'''.rc_validate_pando_evidence_filters <- function(
    padj_threshold, require_padj) {
  if (!is.numeric(padj_threshold) || length(padj_threshold) != 1L ||
      !is.finite(padj_threshold) || padj_threshold < 0 ||
      padj_threshold > 1) {
    stop("`padj_threshold` must be one finite number in [0, 1].",
         call. = FALSE)
  }
  if (!is.logical(require_padj) || length(require_padj) != 1L ||
      is.na(require_padj)) {
    stop("`require_padj` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}'''
)

replace_function(
    'R/pando_grn.R',
    'rc_extract_pando_tf_peak_gene',
    r'''rc_extract_pando_tf_peak_gene <- function(
    grn_object, sample_id, padj_threshold = 0.05,
    require_padj = TRUE) {
  .rc_validate_pando_evidence_filters(
    padj_threshold = padj_threshold, require_padj = require_padj
  )
  if (!is.character(sample_id) || length(sample_id) != 1L ||
      is.na(sample_id) || !nzchar(trimws(sample_id))) {
    stop("`sample_id` must be one non-empty character value.", call. = FALSE)
  }
  if (!requireNamespace("Pando", quietly = TRUE)) {
    stop("Package 'Pando' is required.", call. = FALSE)
  }
  coefs <- as.data.frame(stats::coef(grn_object), stringsAsFactors = FALSE)
  if (!nrow(coefs)) {
    empty <- data.frame(
      sample_id = character(), tf = character(), target = character(),
      region = character(), stringsAsFactors = FALSE
    )
    return(list(all = empty, significant = empty))
  }
  required <- c("tf", "target", "region", "estimate")
  missing <- setdiff(required, colnames(coefs))
  if (length(missing)) {
    stop(
      "Pando coefficient table is missing columns: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  fit <- tryCatch(
    as.data.frame(Pando::gof(grn_object), stringsAsFactors = FALSE),
    error = function(error) data.frame()
  )
  if (nrow(fit) && "target" %in% colnames(fit)) {
    keep_fit <- setdiff(
      colnames(fit),
      intersect(colnames(fit), setdiff(colnames(coefs), "target"))
    )
    coefs <- merge(
      coefs, fit[, keep_fit, drop = FALSE], by = "target",
      all.x = TRUE, sort = FALSE
    )
  }
  coefs$sample_id <- as.character(sample_id)
  coefs$tf <- toupper(as.character(coefs$tf))
  coefs$target <- toupper(as.character(coefs$target))
  coefs$region <- as.character(coefs$region)
  coefs <- coefs[
    , c("sample_id", setdiff(colnames(coefs), "sample_id")), drop = FALSE
  ]
  estimate <- suppressWarnings(as.numeric(coefs$estimate))
  keep <- is.finite(estimate)
  if ("padj" %in% colnames(coefs)) {
    padj <- suppressWarnings(as.numeric(coefs$padj))
    keep <- keep & is.finite(padj) & padj < padj_threshold
  } else if (isTRUE(require_padj)) {
    stop(
      "Pando network does not contain `padj`; use a p-value-producing model.",
      call. = FALSE
    )
  }
  list(all = coefs, significant = coefs[keep, , drop = FALSE])
}'''
)

replace_function(
    'R/standard_pando.R',
    '.rc_filter_standard_pando_edges',
    r'''.rc_filter_standard_pando_edges <- function(table) {
  required <- c("estimate", "padj")
  if (!is.data.frame(table) || !all(required %in% colnames(table))) {
    stop("Standard Pando requires estimate and padj columns.", call. = FALSE)
  }
  estimate <- suppressWarnings(as.numeric(table$estimate))
  padj <- suppressWarnings(as.numeric(table$padj))
  keep <- is.finite(estimate) &
    abs(estimate) > .rc_standard_pando_min_abs_fixed &
    is.finite(padj) & padj < .rc_standard_pando_padj_fixed
  answer <- table[keep, , drop = FALSE]
  attr(answer, "edge_filter") <- list(
    padj = "< 0.05",
    absolute_estimate = paste0(
      "> ", format(.rc_standard_pando_min_abs_fixed, scientific = FALSE)
    )
  )
  answer
}'''
)

standard_path = Path('R/standard_pando.R')
standard_text = standard_path.read_text()
standard_text = replace_once(
    standard_text,
    '''    pando_motif_args = list(), pando_infer_args = list(),
    min_abs_estimate = 0, min_model_rsq = 0.1,
    save_pando_objects = TRUE, parallel = FALSE,
''',
    '''    pando_motif_args = list(), pando_infer_args = list(),
    save_pando_objects = TRUE, parallel = FALSE,
''',
    'standard Pando signature'
)
standard_text = replace_once(
    standard_text,
    '''    extracted <- rc_extract_pando_tf_peak_gene(
      grn_object = grn,
      sample_id = value,
      min_abs_estimate = 0,
      min_model_rsq = min_model_rsq,
      padj_threshold = .rc_standard_pando_padj_fixed,
      require_padj = TRUE
    )
    extracted$significant <- .rc_filter_standard_pando_edges(
      extracted$all,
      min_abs_estimate = min_abs_estimate,
      min_model_rsq = min_model_rsq
    )
''',
    '''    extracted <- rc_extract_pando_tf_peak_gene(
      grn_object = grn,
      sample_id = value,
      padj_threshold = .rc_standard_pando_padj_fixed,
      require_padj = TRUE
    )
    extracted$significant <- .rc_filter_standard_pando_edges(extracted$all)
''',
    'standard Pando extraction filters'
)
standard_text = standard_text.replace(
    '''      predictive_oof_available = FALSE,
      oof_validation_level = "not_applicable_standard_pando",
''',
    ''
)
standard_path.write_text(standard_text)

# Stage 1 accepts only current condition-inference arguments.
step_grn_path = Path('R/step_grn_common_dictionary.R')
step_grn = step_grn_path.read_text()
old_condition_block = '''    standard_only <- intersect(
      names(call_args), c("min_abs_estimate", "min_model_rsq")
    )
    if (length(standard_only)) {
      message(
        "Ignoring standard-Pando-only edge filters in multi-condition mode: ",
        paste(standard_only, collapse = ", "),
        ". The condition penalty is fixed to estimable BH padj < 0.05."
      )
      call_args[standard_only] <- NULL
    }
    retired <- intersect(names(infer_args), c(
      "candidate_screen", "condition_mix", "condition_weight", "alpha",
      "nlambda", "lambda", "lambda_min_ratio", "outer_nfolds",
      "inner_nfolds", "lambda_selection", "scale", "engine_control",
      "comparison_conditions", "active_tol", "max_iter", "tol_objective",
      "tol_coef", "seed", "method", "sample_col", "cv_block_col"
    ))
    if (length(retired)) {
      stop(
        "Retired condition-GRN parameter(s): ",
        paste(retired, collapse = ", "),
        ". Use tf_cor, peak_cor, adjust_method='BH', padj_threshold=0.05, rank_action and min_residual_df.",
        call. = FALSE
      )
    }
'''
new_condition_block = '''    allowed_infer_args <- c(
      "tf_cor", "peak_cor", "adjust_method", "padj_threshold",
      "rank_action", "min_residual_df", "rna_layer", "peak_layer",
      "peak_value_type", "verbose"
    )
    unknown_infer_args <- setdiff(names(infer_args), allowed_infer_args)
    if (length(unknown_infer_args)) {
      stop(
        "Unsupported `pando_infer_args`: ",
        paste(unknown_infer_args, collapse = ", "), call. = FALSE
      )
    }
'''
step_grn = replace_once(
    step_grn, old_condition_block, new_condition_block,
    'Stage 1 old argument cleanup'
)
step_grn_path.write_text(step_grn)

# ---------------------------------------------------------------------------
# Final result exposes current primary and RNA-only control only.
# ---------------------------------------------------------------------------
results_path = Path('R/step_results.R')
results = results_path.read_text()
old_comparisons = '''  common_comparison <- .rc_condition_penalty_route(
    layer2,
    layer2$penalty_common_oof,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  rna_only_comparison <- .rc_condition_penalty_route(
    layer2,
    layer2$penalty_rna_only,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  unique_increment_summary <- .rc_condition_increment_summary(
    layer2,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
'''
new_comparisons = '''  rna_only_comparison <- .rc_condition_penalty_route(
    layer2,
    layer2$penalty_rna_only,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
'''
results = replace_once(
    results, old_comparisons, new_comparisons,
    'result comparison routes'
)
results = replace_once(
    results,
    '''    common_support_component_summary = common_comparison$summary,
    common_support_component_contrast = common_comparison$contrast,
    condition_unique_penalty_increment_summary = unique_increment_summary,
    rna_only_control_summary = rna_only_comparison$summary,
''',
    '''    rna_only_control_summary = rna_only_comparison$summary,
''',
    'result old summaries'
)
results = replace_once(
    results,
    '''      primary_penalty = "condition_full_oof_compatibility_field",
      common_support_role = "compatibility_alias_of_primary",
      condition_unique_role = "zero_compatibility_decomposition",
      nonestimable_edge_policy =
        "coefficient_NA_and_zero_realized_penalty_contribution",
      removed_guardrails = layer2$comparison_contract$removed_guardrails,
''',
    '''      primary_penalty = "penalty",
      rna_only_control = "penalty_rna_only",
      nonestimable_edge_policy =
        "coefficient_NA_and_zero_realized_penalty_contribution",
''',
    'result penalty provenance'
)
results = results.replace(
    'schema_version = "regcompass_regulatory_metabolic_result_v2"',
    'schema_version = "regcompass_regulatory_metabolic_result_v3"'
)
results = results.replace('version = "2.3.0"', 'version = "2.4.0"')
results_path.write_text(results)

# condition_full_contract.R only implemented the removed decomposition summary.
obsolete_source = Path('R/condition_full_contract.R')
if obsolete_source.exists():
    obsolete_source.unlink()

# ---------------------------------------------------------------------------
# DESCRIPTION and current documentation.
# ---------------------------------------------------------------------------
description_path = Path('DESCRIPTION')
description = description_path.read_text()
description = description.replace('Version: 2.3.0', 'Version: 2.4.0')
description = description.replace("    'condition_full_contract.R'\n", '')
description = re.sub(
    r'1667857557/Pando_regcompass@[0-9a-f]{40}',
    f'1667857557/Pando_regcompass@{PANDO_SHA}',
    description
)
description = description.replace(
    'Config/RegCompassR/SupportedSeuratMajors: 4,5',
    'Config/RegCompassR/SupportedSeuratMajors: 4'
)
description = re.sub(
    r'^Config/RegCompassR/SeuratV5MinSignacVersion:.*\n', '',
    description, flags=re.MULTILINE
)
description = description.replace(
    'No pooled-coefficient rescaling, nested condition cross-validation,\n    sparse-group path, structural-zero coefficient substitution or common-\n    support decomposition is used. ',
    'No pooled-coefficient rescaling or alternative condition-effect engine is used. '
)
description_path.write_text(description)

# Remove tests that exclusively validate retired engines or compatibility fields.
obsolete_tests = [
    'tests/testthat/test-v210-oof-biological-contract.R',
    'tests/testthat/test-workflow-hardening-v228.R',
    'tests/testthat/test-detailed-stage-progress.R',
    'tests/testthat/test-layer2-required-mode.R',
    'tests/testthat/test-metadata-contracts.R',
    'tests/testthat/test-stage-io-contracts.R',
    'tests/testthat/test-condition-tutorial-contract.R'
]
for name in obsolete_tests:
    path = Path(name)
    if path.exists():
        path.unlink()

# Update remaining source/test/document references to the current field names.
replacements = {
    'gene_projection_condition_full_oof': 'gene_projection',
    'gene_regulatory_modifier_condition_full_oof': 'gene_regulatory_modifier',
    'gene_support_condition_full_oof': 'gene_support_multiome',
    'reaction_expression_condition_full_oof': 'reaction_expression',
    'reaction_condition_full_support_fraction':
      'reaction_regulatory_support_fraction',
    'penalty_condition_full_oof': 'penalty',
    'score_condition_full_oof_display_only': 'score',
    'score_rna_only_display_only': 'score_rna_only',
    'project_condition_grn_primary_cells': 'project_condition_grn_cells'
}
text_roots = [Path('README.md'), Path('docs'), Path('vignettes'), Path('man'),
              Path('tests/testthat')]
text_files = []
for root in text_roots:
    if root.is_file():
        text_files.append(root)
    elif root.exists():
        text_files.extend(
            path for path in root.rglob('*')
            if path.is_file() and path.suffix.lower() in
            {'.md', '.rmd', '.rd', '.r'}
        )
for path in text_files:
    text = path.read_text()
    for old, new in replacements.items():
        text = text.replace(old, new)
    lines = text.splitlines()
    filtered = []
    for line in lines:
        lower = line.lower()
        if any(token in lower for token in (
            'projection_component', 'comparison_support',
            'regulatory_alpha', 'min_abs_estimate', 'min_model_rsq',
            'penalty_common_oof', 'condition_unique', 'common_oof',
            'nested_cv', 'nested-cv', 'structural_zero',
            'structural-zero', 'predictive_oof', 'oof_validation_level',
            'compatibility_alias', 'compatibility alias'
        )):
            continue
        filtered.append(line)
    path.write_text('\n'.join(filtered).rstrip() + '\n')

# Rewrite current API tests instead of testing deleted fields.
current_test = Path('tests/testthat/test-current-condition-contract.R')
current_test.write_text(r'''test_that("Layer 1 exposes one current projection route", {
  source <- paste(readLines("../../R/layer1_regulatory_support.R"),
                  collapse = "\n")
  expect_match(source, "gene_projection = projection\\$projection")
  expect_match(source, "reaction_expression = reaction_multiome")
  expect_false(grepl("condition_full_oof|common_oof|condition_unique", source))
})

test_that("Layer 2 exposes primary and RNA-only routes", {
  source <- paste(readLines("../../R/step_layer2.R"), collapse = "\n")
  expect_match(source, 'primary = "penalty"', fixed = TRUE)
  expect_match(source, 'rna_control = "penalty_rna_only"', fixed = TRUE)
  expect_false(grepl("penalty_common_oof|penalty_condition_unique", source))
})

test_that("current Pando condition API is used", {
  source <- paste(
    readLines("../../R/condition_grn_contract.R"), collapse = "\n"
  )
  expect_match(source, "Pando::condition_grn_fit(grn_object)", fixed = TRUE)
  expect_match(source, "Pando::project_condition_grn_cells", fixed = TRUE)
  expect_false(grepl("network_name = \\"regcompass_condition_grn\\"", source))
})

test_that("retired user arguments are absent", {
  functions <- list(
    rc_regcompass_step_layer1,
    rc_regcompass_step_grn,
    rc_extract_pando_tf_peak_gene
  )
  arguments <- unique(unlist(lapply(functions, function(fun) names(formals(fun)))))
  expect_false(any(c(
    "projection_component", "comparison_support", "regulatory_alpha",
    "min_abs_estimate", "min_model_rsq"
  ) %in% arguments))
})
''')

# Static final audit.
source_files = sorted(Path('R').glob('*.R'))
source = '\n'.join(path.read_text() for path in source_files)
doc_source = '\n'.join(
    path.read_text() for path in text_files if path.exists()
)
retired_tokens = [
    'condition_full_oof', 'common_oof', 'condition_unique_oof',
    'project_condition_grn_primary_cells',
    'network_name = "regcompass_condition_grn"',
    'min_abs_estimate', 'min_model_rsq', 'structural_zero_probability',
    'nested_cv_used'
]
found = [token for token in retired_tokens if token in source or token in doc_source]
if found:
    raise RuntimeError(f'retired current-contract tokens remain: {found}')
required_tokens = [
    'regcompass_regulatory_layer1_v4',
    'regcompass_regulatory_layer2_v3',
    'gene_projection = projection$projection',
    'reaction_regulatory_support_fraction',
    'primary = "penalty"',
    'rna_control = "penalty_rna_only"',
    PANDO_SHA
]
missing = [token for token in required_tokens
           if token not in source and token not in description]
if missing:
    raise RuntimeError(f'current contract tokens missing: {missing}')
