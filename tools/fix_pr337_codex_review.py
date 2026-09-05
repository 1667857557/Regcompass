from pathlib import Path


def replace_once(path, old, new):
    p = Path(path)
    text = p.read_text()
    count = text.count(old)
    if count != 1:
        raise SystemExit(f"{path}: expected one match, found {count}: {old[:120]!r}")
    p.write_text(text.replace(old, new, 1))


# Package version for the review-fix revision.
replace_once(
    "DESCRIPTION",
    "Version: 2.4.24\n",
    "Version: 2.4.25\n",
)

# 1) Keep paired Layer 2 progress monotonic: phase 6 must precede phase 8.
replace_once(
    "R/step_layer2.R",
    '''  if (!is.null(args$control_layer1)) {\n    .rc_layer2_overall_event(\n      "rna_control_complete", 8L,\n      "RNA-only Step 2 completed in the same target dispatch with shared Vmax"\n    )\n  }\n  .rc_layer2_overall_event(\n    "primary_engine_complete", 6L,\n    "primary structural models and directional scores assembled"\n  )\n''',
    '''  .rc_layer2_overall_event(\n    "primary_engine_complete", 6L,\n    "primary structural models and directional scores assembled"\n  )\n  if (!is.null(args$control_layer1)) {\n    .rc_layer2_overall_event(\n      "rna_control_complete", 8L,\n      "RNA-only Step 2 completed in the same target dispatch with shared Vmax"\n    )\n  }\n''',
)

# 2) The two objective streams are independent LP solves but share one
# persistent target solver. Keep the legacy field but make its meaning correct.
replace_once(
    "R/step_layer2.R",
    '''    paired_step2_dispatch = TRUE,\n    independent_solver_streams = TRUE,\n    exact_identical_model_objective_reuse = TRUE,\n''',
    '''    paired_step2_dispatch = TRUE,\n    independent_solver_streams = FALSE,\n    independent_lp_solves_on_shared_target_engine = TRUE,\n    exact_identical_model_objective_reuse = TRUE,\n''',
)

# 3) Preserve the structural model's target status. Solver feasibility/status
# remains separately reported and must not overwrite reconstruction provenance.
replace_once(
    "R/celltype_microcompass_reaction_parallel.R",
    '''    diagnostics <- data.frame(\n      row_id = rep(row_id, length(units)),\n''',
    '''    target_status <- if (!is.null(model$target_status)) {\n      rep(as.character(model$target_status), length(units))\n    } else {\n      ifelse(primary$feasible, "ok", "structurally_infeasible")\n    }\n    diagnostics <- data.frame(\n      row_id = rep(row_id, length(units)),\n''',
)
replace_once(
    "R/celltype_microcompass_reaction_parallel.R",
    '''      target_status = ifelse(\n        primary$feasible, "ok", "structurally_infeasible"\n      ),\n''',
    '''      target_status = target_status,\n''',
)
replace_once(
    "R/celltype_microcompass_reaction_parallel.R",
    '''      control_diagnostics, prepared, primary, control, primary_metrics,\n      control_metrics\n''',
    '''      control_diagnostics, prepared, primary, control, primary_metrics,\n      control_metrics, target_status\n''',
)

# 4) Reuse flags must be keyed by target row because checkpoint order is
# model/batch order, not necessarily the original penalty row order. Apply the
# same protection to both full-GEM and cell-type engines.
for path, label in [
    ("R/microcompass_engine.R", "full-GEM"),
    ("R/celltype_microcompass_reaction_parallel.R", "cell-type"),
]:
    replace_once(
        path,
        '''  control_reused <- logical(length(checkpoint_files))\n  observed_rows <- character(length(checkpoint_files))\n''',
        '''  control_reused <- setNames(logical(length(row_ids)), row_ids)\n  observed_rows <- character(length(checkpoint_files))\n''',
    )
    replace_once(
        path,
        '''      control_reused[[i]] <- isTRUE(all(reuse_mask))\n''',
        '''      control_reused[[row_id]] <- isTRUE(all(reuse_mask))\n''',
    )
    replace_once(
        path,
        '''    rna_control_model_identical_reuse = control_reused,\n''',
        '''    rna_control_model_identical_reuse = control_reused[row_ids],\n''',
    )

# Update misleading per-engine prose in the full-GEM result metadata.
replace_once(
    "R/microcompass_engine.R",
    '''      step2_solver_reuse = paste(\n        "one prepared target template with independent persistent HiGHS",\n        "streams for primary and RNA-only objectives across all metacells"\n      ),\n''',
    '''      step2_solver_reuse = paste(\n        "one prepared target template and one persistent HiGHS engine per",\n        "target; primary and RNA-only remain independent LP solves when their",\n        "full model-wide objective vectors differ"\n      ),\n''',
)

# Extend the permanent paired-control regression.
replace_once(
    "tests/layer2-paired-control-check.R",
    '''mk_payload <- function(file, primary, control = NULL, full = TRUE,\n                       direction = "forward") {\n''',
    '''mk_payload <- function(file, primary, control = NULL, full = TRUE,\n                       direction = "forward", model_status = "ok") {\n''',
)
replace_once(
    "tests/layer2-paired-control-check.R",
    '''      S = S, lb = lb, ub = ub, target_status = "ok",\n''',
    '''      S = S, lb = lb, ub = ub, target_status = model_status,\n''',
)
replace_once(
    "tests/layer2-paired-control-check.R",
    '''run_worker <- function(worker, primary, control = NULL, full = TRUE,\n                       direction = "forward") {\n''',
    '''run_worker <- function(worker, primary, control = NULL, full = TRUE,\n                       direction = "forward", model_status = "ok") {\n''',
)
replace_once(
    "tests/layer2-paired-control-check.R",
    '''  mk_payload(payload, primary, control, full, direction)\n''',
    '''  mk_payload(payload, primary, control, full, direction, model_status)\n''',
)

insert_before = '''for (direction in c("forward", "reverse")) {\n  check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE, direction)\n'''
insert_text = '''# Structural target provenance must not be rewritten from Step 2 solver\n# feasibility. A deliberately distinctive model status must survive unchanged.\nstatus_check <- run_worker(\n  .rc_step2_reaction_batch_worker, primary_penalty,\n  full = FALSE, model_status = "structural_status_preserved"\n)\nstopifnot(all(\n  status_check$diagnostics$target_status == "structural_status_preserved"\n))\n\nfor (direction in c("forward", "reverse")) {\n  check_worker(.rc_full_gem_step2_reaction_batch_worker, TRUE, direction)\n'''
replace_once("tests/layer2-paired-control-check.R", insert_before, insert_text)

replace_once(
    "tests/layer2-paired-control-check.R",
    '''  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),\n''',
    '''  grepl("paired_step2_dispatch = TRUE", stage, fixed = TRUE),\n  grepl("independent_solver_streams = FALSE", stage, fixed = TRUE),\n  grepl("independent_lp_solves_on_shared_target_engine = TRUE", stage,\n        fixed = TRUE),\n''',
)
replace_once(
    "tests/layer2-paired-control-check.R",
    '''  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)\n)\n''',
    '''  grepl("answer$comparison_paths <- NULL", stage, fixed = TRUE)\n)\n\n# The paired progress wrapper must never emit phase 8 before phase 6.\nphase6 <- regexpr('"primary_engine_complete", 6L', stage, fixed = TRUE)[[1L]]\nphase8 <- regexpr('"rna_control_complete", 8L', stage, fixed = TRUE)[[1L]]\nstopifnot(phase6 > 0L, phase8 > 0L, phase6 < phase8)\n\n# Reuse summaries are row-keyed in both engines, so model-scoped checkpoint\n# ordering cannot scramble target-level diagnostics.\nfull_text <- paste(readLines("R/microcompass_engine.R", warn = FALSE),\n                   collapse = "\\n")\ncell_text <- paste(readLines(\n  "R/celltype_microcompass_reaction_parallel.R", warn = FALSE\n), collapse = "\\n")\nfor (text in list(full_text, cell_text)) {\n  stopifnot(\n    grepl("control_reused <- setNames(logical(length(row_ids)), row_ids)",\n          text, fixed = TRUE),\n    grepl("control_reused[[row_id]] <- isTRUE(all(reuse_mask))",\n          text, fixed = TRUE),\n    grepl("rna_control_model_identical_reuse = control_reused[row_ids]",\n          text, fixed = TRUE)\n  )\n}\n''',
)

print("Applied PR #337 Codex review fixes")
