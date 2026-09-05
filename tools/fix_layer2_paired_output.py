from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def patch(path, old, new, label):
    p = ROOT / path
    text = p.read_text()
    if old not in text:
        raise RuntimeError(f"{label}: pattern not found in {path}")
    p.write_text(text.replace(old, new, 1))

# The meta-module RNA-control direction contract requires cell_type provenance.
patch(
    "R/celltype_microcompass_reaction_parallel.R",
    '''      data.frame(\n        row_id = rep(row_id, length(units)),\n        unit_id = units,\n        reaction_id = rep(entry$reaction_id, length(units)),\n        target_direction = rep(entry$target_direction, length(units)),''',
    '''      data.frame(\n        row_id = rep(row_id, length(units)),\n        unit_id = units,\n        cell_type = rep(entry$cell_type, length(units)),\n        reaction_id = rep(entry$reaction_id, length(units)),\n        target_direction = rep(entry$target_direction, length(units)),''',
    "cell-type control diagnostics"
)

# Do not retain a near-complete second microCOMPASS result merely to validate
# the RNA-control contract. The persisted storage policy already discards it;
# validate with the minimal fields actually read by the contract instead.
patch(
    "R/step_layer2.R",
    '''  rna_only <- answer\n  rna_only$penalty <- answer$penalty_rna_only\n  rna_only$score <- answer$score_rna_only\n  rna_only$feasible <- answer$feasible_rna_only\n  rna_only$evaluated <- answer$evaluated_rna_only\n  rna_only$lp_diagnostics <- answer$lp_diagnostics_rna_only\n  rna_only$step2_engine_metrics <- answer$step2_engine_metrics_rna_only\n  rna_only$shared_model_cache <- NULL\n  .rc_assert_layer2_shared_contract(answer, rna_only, "RNA-only control")\n\n  answer$schema_version''',
    '''  rna_only_contract <- list(\n    penalty = answer$penalty_rna_only,\n    vmax = answer$vmax,\n    feasible = answer$feasible_rna_only,\n    evaluated = answer$evaluated_rna_only,\n    structural_model_contract = answer$structural_model_contract,\n    lp_diagnostics = answer$lp_diagnostics_rna_only,\n    model_mode = answer$model_mode,\n    medium_scenarios = answer$medium_scenarios,\n    unit_meta = answer$unit_meta\n  )\n  .rc_assert_layer2_shared_contract(\n    answer, rna_only_contract, "RNA-only control"\n  )\n  rm(rna_only_contract)\n\n  answer$schema_version''',
    "lightweight RNA-control contract"
)
patch(
    "R/step_layer2.R",
    '''  answer$comparison_paths <- list(rna_only = rna_only)''',
    '''  answer$comparison_paths <- NULL''',
    "drop transient duplicate comparison result"
)

# The standalone numerical regression deliberately sources only the hot LP files,
# matching the existing native regression. Provide the same bound-alignment helper
# that package collation supplies in an installed namespace.
patch(
    "tests/layer2-paired-control-check.R",
    '''source("R/00_utils.R", local = FALSE)\nsource("R/lp_solver.R", local = FALSE)''',
    '''source("R/00_utils.R", local = FALSE)\n\nrc_align_bound <- function(x, reactions, default, name) {\n  if (is.null(x)) {\n    value <- rep(default, length(reactions))\n    names(value) <- reactions\n    return(value)\n  }\n  if (!is.null(names(x))) {\n    missing <- setdiff(reactions, names(x))\n    if (length(missing)) {\n      stop(name, " is missing reactions: ", paste(missing, collapse = ", "))\n    }\n    x <- x[reactions]\n  }\n  x <- as.numeric(x)\n  if (length(x) != length(reactions) || anyNA(x)) {\n    stop(name, " does not align with reactions")\n  }\n  names(x) <- reactions\n  x\n}\n\nsource("R/layer2_parallel_runtime.R", local = FALSE)\nsource("R/lp_solver.R", local = FALSE)''',
    "standalone LP helpers"
)
patch(
    "tests/layer2-paired-control-check.R",
    '''  !grepl("rna_only <- run_control(", stage, fixed = TRUE)\n)''',
    '''  !grepl("rna_only <- run_control(", stage, fixed = TRUE),\n  !grepl("rna_only <- answer", stage, fixed = TRUE),\n  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)\n)''',
    "assert lightweight stage result"
)

# Keep the structural test in sync with the memory-safe compatibility path.
patch(
    "tests/testthat/test-layer2-paired-control.R",
    '''  expect_false(grepl("rna_only <- run_control(", stage, fixed = TRUE))''',
    '''  expect_false(grepl("rna_only <- run_control(", stage, fixed = TRUE))\n  expect_false(grepl("rna_only <- answer", stage, fixed = TRUE))\n  expect_match(stage, "answer$comparison_paths <- NULL", fixed = TRUE)''',
    "test lightweight stage result"
)

print("Paired Layer 2 output patch tightened")
