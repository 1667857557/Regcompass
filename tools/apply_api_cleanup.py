from __future__ import annotations

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
    target = ROOT / path
    target.write_text(text, encoding="utf-8")


def function_matches(text: str, name: str) -> list[re.Match[str]]:
    pattern = re.compile(
        rf"(?m)^(?P<indent>[ \t]*){re.escape(name)}[ \t]*<-[ \t]*function[ \t]*\("
    )
    return list(pattern.finditer(text))


def find_function_span(text: str, name: str, occurrence: int = 1) -> tuple[int, int]:
    matches = function_matches(text, name)
    if len(matches) < occurrence:
        raise RuntimeError(
            f"Expected function {name!r} occurrence {occurrence}; found {len(matches)}"
        )
    start = matches[occurrence - 1].start()
    pos = matches[occurrence - 1].end()
    quote: str | None = None
    escaped = False
    comment = False
    body_start = None
    while pos < len(text):
        char = text[pos]
        if comment:
            if char == "\n":
                comment = False
            pos += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\" and quote != "`":
                escaped = True
            elif char == quote:
                quote = None
            pos += 1
            continue
        if char == "#":
            comment = True
        elif char in ('"', "'", "`"):
            quote = char
        elif char == "{":
            body_start = pos
            break
        pos += 1
    if body_start is None:
        raise RuntimeError(f"Could not locate body for {name!r}")

    depth = 0
    quote = None
    escaped = False
    comment = False
    pos = body_start
    while pos < len(text):
        char = text[pos]
        if comment:
            if char == "\n":
                comment = False
            pos += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\" and quote != "`":
                escaped = True
            elif char == quote:
                quote = None
            pos += 1
            continue
        if char == "#":
            comment = True
        elif char in ('"', "'", "`"):
            quote = char
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                end = pos + 1
                while end < len(text) and text[end] in " \t":
                    end += 1
                if end < len(text) and text[end] == "\n":
                    end += 1
                while end < len(text) and text[end] == "\n":
                    end += 1
                return start, end
        pos += 1
    raise RuntimeError(f"Unclosed body for {name!r}")


def remove_function(path: str, name: str, occurrence: int = 1) -> None:
    text = read(path)
    start, end = find_function_span(text, name, occurrence)
    write(path, text[:start] + text[end:])


def rename_function(path: str, old: str, new: str, occurrence: int = 1) -> None:
    text = read(path)
    matches = function_matches(text, old)
    if len(matches) < occurrence:
        raise RuntimeError(
            f"Expected function {old!r} occurrence {occurrence}; found {len(matches)}"
        )
    match = matches[occurrence - 1]
    name_start = match.start() + len(match.group("indent"))
    write(path, text[:name_start] + new + text[name_start + len(old):])


def remove_alias(path: str, alias: str, target: str) -> None:
    text = read(path)
    pattern = re.compile(
        rf"(?m)^[ \t]*{re.escape(alias)}[ \t]*<-[ \t]*(?:\n[ \t]*)?"
        rf"{re.escape(target)}[ \t]*\n+"
    )
    updated, count = pattern.subn("", text, count=1)
    if count != 1:
        raise RuntimeError(
            f"Expected one alias {alias!r} <- {target!r} in {path}; found {count}"
        )
    write(path, updated)


def replace(path: str, old: str, new: str, count: int | None = None) -> None:
    text = read(path)
    observed = text.count(old)
    if observed == 0:
        raise RuntimeError(f"Did not find {old!r} in {path}")
    if count is not None and observed != count:
        raise RuntimeError(
            f"Expected {count} occurrences of {old!r} in {path}; found {observed}"
        )
    write(path, text.replace(old, new))


def remove_line(path: str, exact: str) -> None:
    text = read(path)
    line = exact + "\n"
    if text.count(line) != 1:
        raise RuntimeError(f"Expected exactly one line {exact!r} in {path}")
    write(path, text.replace(line, ""))


def stage1_cleanup() -> None:
    rename_function(
        "R/condition_grn_contract.R",
        ".rc_require_pando_condition_grn_fit",
        ".rc_require_pando_condition_grn_fit_schema",
    )
    rename_function(
        "R/condition_grn_contract.R",
        ".rc_extract_condition_grn_contract",
        ".rc_extract_condition_grn_contract_core",
    )
    remove_alias(
        "R/condition_grn_contract_compat.R",
        ".rc_require_pando_condition_grn_fit_base",
        ".rc_require_pando_condition_grn_fit",
    )
    remove_alias(
        "R/condition_grn_contract_compat.R",
        ".rc_extract_condition_grn_contract_impl",
        ".rc_extract_condition_grn_contract",
    )
    replace(
        "R/condition_grn_contract_compat.R",
        ".rc_require_pando_condition_grn_fit_base",
        ".rc_require_pando_condition_grn_fit_schema",
    )
    replace(
        "R/condition_grn_contract_compat.R",
        ".rc_extract_condition_grn_contract_impl",
        ".rc_extract_condition_grn_contract_core",
    )

    rename_function(
        "R/mixed_pando.R",
        ".rc_merge_pando_results",
        ".rc_merge_pando_results_core",
    )
    remove_function("R/mixed_pando.R", ".rc_standard_pando_runtime_args")
    remove_function("R/mixed_pando.R", ".rc_fit_pando_by_celltype_route")
    remove_function("R/mixed_pando.R", ".rc_overlay_projection")

    remove_alias(
        "R/mixed_pando_contract.R",
        ".rc_merge_pando_results_contract_base",
        ".rc_merge_pando_results",
    )
    rename_function(
        "R/mixed_pando_contract.R",
        ".rc_merge_pando_results",
        ".rc_merge_pando_results_validated",
    )
    replace(
        "R/mixed_pando_contract.R",
        ".rc_merge_pando_results_contract_base",
        ".rc_merge_pando_results_core",
    )

    rename_function(
        "R/grn_route_parallel.R",
        ".rc_route_pando_infer_args",
        ".rc_route_pando_infer_args_core",
    )
    remove_function("R/grn_route_parallel.R", ".rc_merge_pando_grn_data")
    remove_function("R/grn_route_parallel.R", ".rc_merge_condition_job_results")

    remove_alias(
        "R/grn_integration_contract.R",
        ".rc_route_pando_infer_args_integration_base",
        ".rc_route_pando_infer_args",
    )
    replace(
        "R/grn_integration_contract.R",
        ".rc_route_pando_infer_args_integration_base",
        ".rc_route_pando_infer_args_core",
    )

    rename_function(
        "R/grn_parallel_output_contract.R",
        ".rc_merge_condition_job_results",
        ".rc_merge_condition_job_results_core",
    )
    remove_alias(
        "R/grn_parallel_output_contract.R",
        ".rc_merge_pando_results_base_parallel_contract",
        ".rc_merge_pando_results",
    )
    rename_function(
        "R/grn_parallel_output_contract.R",
        ".rc_merge_pando_results",
        ".rc_merge_pando_results_with_parallel_objects",
    )
    replace(
        "R/grn_parallel_output_contract.R",
        ".rc_merge_pando_results_base_parallel_contract",
        ".rc_merge_pando_results_validated",
    )

    remove_alias(
        "R/grn_parallel_merge_contract.R",
        ".rc_merge_condition_job_results_contract_base",
        ".rc_merge_condition_job_results",
    )
    replace(
        "R/grn_parallel_merge_contract.R",
        ".rc_merge_condition_job_results_contract_base",
        ".rc_merge_condition_job_results_core",
    )

    remove_alias(
        "R/condition_penalty_contract.R",
        ".rc_merge_pando_results_penalty_base",
        ".rc_merge_pando_results",
    )
    replace(
        "R/condition_penalty_contract.R",
        ".rc_merge_pando_results_penalty_base",
        ".rc_merge_pando_results_with_parallel_objects",
    )

    remove_line("R/standard_pando.R", ".rc_standard_pando_min_abs_fixed <- 0.01")
    remove_function("R/standard_pando.R", ".rc_filter_standard_pando_edges")
    replace(
        "R/standard_pando.R",
        ".rc_standard_pando_min_abs_fixed",
        ".RC_PANDO_PENALTY_ESTIMATE_THRESHOLD",
    )


def remove_confirmed_dead_functions() -> None:
    dead = {
        "R/00_utils.R": ["rc_drop_na_grouping", ".rc_clamp01"],
        "R/execution_monitor.R": [".rc_write_execution_timing"],
        "R/gpr_capacity.R": ["rc_robust_z", ".rc_absolute_activity_score"],
        "R/layer2_corda_direction_contract.R": [".rc_corda_results_table"],
        "R/layer2_corda_evidence.R": [".rc_corda_integer"],
        "R/layer2_corda_lp.R": [
            ".rc_corda_block_reactions",
            ".rc_corda_task_bpparam",
            ".rc_corda_worker_count",
        ],
        "R/layer2_corda_target_contract.R": [
            ".rc_corda_apply_target_epsilon"
        ],
        "R/metacell_fragments.R": [
            ".rc_register_signac_fragments",
            ".rc_restore_metacell_metadata",
        ],
        "R/microcompass_engine.R": ["rc_summarize_microcompass"],
    }
    for path, names in dead.items():
        for name in names:
            remove_function(path, name)


def update_readme() -> None:
    path = "README.md"
    text = read(path)
    heading = "## Simple workflow\n"
    if heading not in text:
        raise RuntimeError("README simple workflow heading not found")
    diagram = """## Workflow\n\n```mermaid\nflowchart TD\n  A[Paired single-cell RNA + ATAC Seurat object] --> S1\n  G[Species GEM] --> S3\n  M[Shared biological or custom medium] --> S5\n\n  S1[Stage 1: cell-type Pando GRN routing] --> S2[Stage 2: condition-pure multimodal metacells]\n  S2 --> S3[Stage 3: cell-type reaction meta-modules]\n  S1 --> S4[Stage 4: RNA + regulatory reaction support]\n  S2 --> S4\n  S3 --> S4\n\n  S3 --> S5[Stage 5 default: CORDA2 cell-type x medium structural model]\n  S4 --> P[Primary multiome COMPASS-like directional penalty]\n  S4 --> R[RNA-only control penalty]\n  S5 --> P\n  S5 --> R\n\n  P --> S6[Stage 6: rankings, annotations and condition contrasts]\n  R --> S6\n\n  F[Explicit supplement: FASTCORE] -. replaces CORDA2 .-> S5\n  H[Explicit supplement: complete full GEM] -. skips reconstruction .-> P\n```\n\nThe primary multiome score and RNA-only control reuse the same completed structural model, bounds, medium and target directions. Therefore their difference isolates the regulatory contribution rather than a change in network structure.\n\n"""
    if "## Workflow\n" not in text:
        text = text.replace(heading, diagram + heading, 1)
    write(path, text)


def update_description() -> None:
    path = "DESCRIPTION"
    text = read(path)
    old = (
        "    FASTCORE is the default compact add-only completion method; original MATLAB\n"
        "    CORDA2 and a COMPASS-style full-GEM route are explicit alternatives."
    )
    new = (
        "    Original MATLAB CORDA2 is the default context-specific completion method;\n"
        "    FASTCORE and a COMPASS-style full-GEM route are explicit supplements."
    )
    if old not in text:
        raise RuntimeError("DESCRIPTION Layer 2 default text not found")
    write(path, text.replace(old, new, 1))


def main() -> None:
    stage1_cleanup()
    remove_confirmed_dead_functions()
    update_readme()
    update_description()


if __name__ == "__main__":
    main()
