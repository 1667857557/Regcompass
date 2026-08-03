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


layer1_validator = r'''.rc_validate_layer1_stage <- function(
    layer1, workflow_params = NULL, gem = NULL, argument = "layer1") {
  .rc_require_stage_class(
    layer1, "regcompass_layer1_step", argument, "rc_regcompass_step_layer1"
  )
  ids <- .rc_layer1_unit_ids(layer1)
  mode <- layer1$analysis_mode
  valid_modes <- c("condition_grn", "standard_pando", "mixed_pando")
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% valid_modes) {
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
    stop("Layer 1 regulatory modifier is inconsistent with its inputs.",
         call. = FALSE)
  }
  value <- layer1$reaction_expression_rna_only
  if (!is.numeric(value) || is.null(dim(value)) ||
      !identical(dimnames(value), dimnames(layer1$reaction_expression))) {
    stop("Layer 1 RNA-only reaction matrix is misaligned.", call. = FALSE)
  }
  if (!is.numeric(layer1$reaction_regulatory_support_fraction) ||
      !identical(
        dimnames(layer1$reaction_regulatory_support_fraction),
        dimnames(layer1$reaction_expression)
      )) {
    stop("Layer 1 reaction regulatory support is misaligned.",
         call. = FALSE)
  }

  provenance <- layer1$projection_provenance
  required_provenance <- c(
    "analysis_mode", "pando_schema", "projection_origin",
    "projection_used_for_penalty", "projection_name",
    "condition_coefficients_calculated", "supercell_membership",
    "unavailable_target_policy", "nonestimable_edge_policy",
    "cell_type_analysis_mode"
  )
  routing <- provenance$cell_type_analysis_mode
  if (!identical(layer1$schema_version, "regcompass_regulatory_layer1_v4") ||
      !is.list(provenance) ||
      !all(required_provenance %in% names(provenance)) ||
      !identical(provenance$analysis_mode, mode) ||
      !isTRUE(provenance$projection_used_for_penalty) ||
      !identical(
        provenance$supercell_membership,
        "membership_table(cell_id, metacell_id)"
      ) ||
      !is.data.frame(routing) ||
      !all(c("cell_type", "analysis_mode") %in% colnames(routing)) ||
      !nrow(routing) ||
      anyNA(routing$cell_type) || anyNA(routing$analysis_mode) ||
      any(!routing$analysis_mode %in% c("condition_grn", "standard_pando")) ||
      any(!nzchar(as.character(provenance$pando_schema))) ||
      any(!nzchar(as.character(provenance$projection_origin))) ||
      any(!nzchar(as.character(provenance$projection_name))) ||
      any(!nzchar(as.character(provenance$nonestimable_edge_policy)))) {
    stop("Layer 1 regulatory projection provenance is incompatible.",
         call. = FALSE)
  }
  observed_modes <- sort(unique(as.character(routing$analysis_mode)))
  expected_modes <- switch(
    mode,
    condition_grn = "condition_grn",
    standard_pando = "standard_pando",
    mixed_pando = sort(c("condition_grn", "standard_pando"))
  )
  if (!identical(observed_modes, expected_modes)) {
    stop("Layer 1 cell-type routing does not match its analysis mode.",
         call. = FALSE)
  }
  expected_coefficients <- "condition_grn" %in% observed_modes
  if (!identical(
        provenance$condition_coefficients_calculated,
        expected_coefficients
      )) {
    stop("Layer 1 condition-coefficient provenance is inconsistent.",
         call. = FALSE)
  }
  if (!is.null(workflow_params)) {
    .rc_require_workflow_params(layer1, workflow_params, argument)
  }
  if (!is.null(gem)) .rc_require_stage_gem(layer1, gem, argument)
  invisible(ids)
}'''

layer2_validator = r'''.rc_validate_layer2_stage <- function(
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
      stop("`required_mode` must be one of: ",
           paste(valid_modes, collapse = ", "), ".", call. = FALSE)
    }
    if (!identical(layer2$model_mode, required_mode)) {
      stop("Layer 2 model mode is `", layer2$model_mode,
           "`; `", required_mode, "` is required.", call. = FALSE)
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

replace_function("R/stage_contracts.R",
                 ".rc_validate_layer1_stage", layer1_validator)
replace_function("R/stage_contracts.R",
                 ".rc_validate_layer2_stage", layer2_validator)

step = Path("R/step_grn_common_dictionary.R")
text = step.read_text()
text = text.replace(
'''      "rank_action", "min_residual_df", "rna_layer", "peak_layer",
      "peak_value_type", "verbose"
''',
'''      "rank_action", "min_residual_df", "rna_layer", "peak_layer",
      "peak_value_type"
''', 1)
text = text.replace(
'''      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L,
      verbose = .rc_progress_enabled(progress)
''',
'''      padj_threshold = 0.05, rank_action = "mark",
      min_residual_df = 1L
''', 1)
step.write_text(text)

for path in Path("R").glob("*.R"):
    text = path.read_text()
    text = text.replace("outer-heldout", "condition-specific")
    text = text.replace("OOF", "condition")
    text = text.replace("oof", "condition")
    path.write_text(text)

source = "\n".join(p.read_text() for p in Path("R").glob("*.R"))
for token in (
    "condition_full_oof", "common_oof", "condition_unique_oof",
    'network_name = "regcompass_condition_grn"',
    "min_model_rsq", "min_abs_estimate",
    "predictive_oof_available", "oof_validation_level"
):
    if token in source:
        raise RuntimeError(f"retired token remains in R source: {token}")
