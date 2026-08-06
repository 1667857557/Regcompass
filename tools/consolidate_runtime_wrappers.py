from __future__ import annotations

from pathlib import Path
import re

import apply_api_cleanup as util

ROOT = Path(__file__).resolve().parents[1]
POOL = "R/layer2_corda_pool_lifecycle.R"


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def write(path: str, text: str) -> None:
    (ROOT / path).write_text(text.rstrip() + "\n", encoding="utf-8")


def function_text(path: str, name: str) -> str:
    text = read(path)
    start, end = util.find_function_span(text, name)
    return text[start:end].rstrip()


def remove_function(path: str, name: str) -> None:
    text = read(path)
    start, end = util.find_function_span(text, name)
    write(path, text[:start] + text[end:])


def rename_function(path: str, old: str, new: str) -> None:
    text = read(path)
    matches = util.function_matches(text, old)
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one definition of {old} in {path}; found {len(matches)}"
        )
    match = matches[0]
    start = match.start() + len(match.group("indent"))
    write(path, text[:start] + new + text[start + len(old):])


def remove_capture(name: str) -> None:
    text = read(POOL)
    pattern = re.compile(
        rf"(?m)^[ \t]*[.]rc_layer2_progress_capture\(\"{re.escape(name)}\"\)"
        rf"[ \t]*\n"
    )
    updated, count = pattern.subn("", text, count=1)
    if count != 1:
        raise RuntimeError(f"Expected one progress capture for {name}")
    write(POOL, updated)


def append_wrapper(target_file: str, wrapper: str, label: str) -> None:
    text = read(target_file).rstrip()
    addition = (
        "\n\n# Progress-aware entry point; the algorithm remains in the core above.\n"
        + wrapper.rstrip()
        + "\n"
    )
    if label in text:
        raise RuntimeError(f"Wrapper {label} already exists in {target_file}")
    write(target_file, text + addition)


def move_pool_wrapper(name: str, target_file: str, core_name: str) -> None:
    wrapper = function_text(POOL, name)
    expected = f".rc_layer2_progress_original${name}"
    if expected not in wrapper:
        raise RuntimeError(f"Pool wrapper {name} does not call its captured implementation")
    wrapper = wrapper.replace(expected, core_name)
    rename_function(target_file, name, core_name)
    remove_capture(name)
    remove_function(POOL, name)
    append_wrapper(target_file, wrapper, f"{name} <- function")


def consolidate_pool_wrappers() -> None:
    mapping = [
        (".rc_build_celltype_medium_union_gem_cache", "R/celltype_union_gem_cache.R", ".rc_build_celltype_medium_union_gem_cache_core"),
        (".rc_build_microcompass_vmax_cache", "R/microcompass_vmax_cache.R", ".rc_build_microcompass_vmax_cache_core"),
        (".rc_complete_celltype_medium_corda_gem", "R/layer2_corda_model.R", ".rc_complete_celltype_medium_corda_gem_core"),
        (".rc_complete_celltype_medium_union_gem", "R/fastcore.R", ".rc_complete_celltype_medium_union_gem_core"),
        (".rc_corda_build_three_stage", "R/layer2_corda2_algorithm_build.R", ".rc_corda_build_three_stage_core"),
        (".rc_corda_classify_reactions", "R/layer2_corda_evidence.R", ".rc_corda_classify_reactions_core"),
        (".rc_corda_core_closure", "R/layer2_corda_model.R", ".rc_corda_core_closure_core"),
        (".rc_corda_parent", "R/layer2_corda_parent_contract.R", ".rc_corda_parent_core"),
        (".rc_corda2_dependency_assessment", "R/layer2_corda2_algorithm.R", ".rc_corda2_dependency_assessment_core"),
        (".rc_corda2_maximize_target", "R/layer2_corda_paper_contract.R", ".rc_corda2_maximize_target_core"),
        (".rc_directional_feasibility", "R/fastcore.R", ".rc_directional_feasibility_core"),
        (".rc_fastcore_complete_direction", "R/fastcore.R", ".rc_fastcore_complete_direction_core"),
        (".rc_fastcore_parent", "R/fastcore.R", ".rc_fastcore_parent_core"),
        (".rc_layer2_comparison_table", "R/step_layer2.R", ".rc_layer2_comparison_table_core"),
        (".rc_run_microcompass_engine_monitored", "R/step_layer2.R", ".rc_run_microcompass_engine_monitored_core"),
        (".rc_run_microcompass_monitored", "R/step_layer2.R", ".rc_run_microcompass_monitored_core"),
        (".rc_run_shared_full_gem_engine", "R/microcompass_engine.R", ".rc_run_shared_full_gem_engine_core"),
        (".rc_step_monitor_fail", "R/execution_monitor.R", ".rc_step_monitor_fail_core"),
        (".rc_step_monitor_finish", "R/execution_monitor.R", ".rc_step_monitor_finish_core"),
        (".rc_step_monitor_start", "R/execution_monitor.R", ".rc_step_monitor_start_core"),
        (".rc_validate_layer2_stage", "R/stage_contracts.R", ".rc_validate_layer2_stage_core"),
        ("rc_build_full_gem", "R/full_gem.R", ".rc_build_full_gem_core"),
        ("rc_build_full_gem_cache", "R/full_gem.R", ".rc_build_full_gem_cache_core"),
    ]
    for name, target, core in mapping:
        move_pool_wrapper(name, target, core)

    # The reaction-parallel implementation supersedes the older cell-type engine.
    remove_function("R/celltype_microcompass_engine.R", ".rc_run_celltype_microcompass_engine")
    rename_function(
        "R/celltype_microcompass_reaction_parallel.R",
        ".rc_run_celltype_microcompass_engine",
        ".rc_run_celltype_microcompass_engine_reaction_core",
    )
    wrapper = function_text(POOL, ".rc_run_celltype_microcompass_engine")
    expected = ".rc_layer2_progress_original$.rc_run_celltype_microcompass_engine"
    if expected not in wrapper:
        raise RuntimeError("Reaction scoring wrapper lost its captured implementation")
    wrapper = wrapper.replace(
        expected, ".rc_run_celltype_microcompass_engine_reaction_core"
    )
    remove_capture(".rc_run_celltype_microcompass_engine")
    remove_function(POOL, ".rc_run_celltype_microcompass_engine")
    append_wrapper(
        "R/celltype_microcompass_reaction_parallel.R",
        wrapper,
        ".rc_run_celltype_microcompass_engine <- function",
    )

    pool = read(POOL)
    pool = re.sub(
        r"(?m)^[.]rc_layer2_progress_original <- new[.]env\(parent = emptyenv\(\)\)\n",
        "",
        pool,
        count=1,
    )
    write(POOL, pool)
    remove_function(POOL, ".rc_layer2_progress_capture")
    pool = read(POOL)
    if ".rc_layer2_progress_original" in pool or ".rc_layer2_progress_capture" in pool:
        raise RuntimeError("Runtime capture infrastructure remains after consolidation")
    pool = pool.replace(
        "# This section is loaded after every Layer 2 implementation so the wrappers\n# below observe the final runtime functions, including the reaction-parallel\n# scoring engine and the original CORDA2 implementation.\n",
        "# Canonical Layer 2 functions call these shared progress helpers directly.\n# No runtime function capture or same-name redefinition is used.\n",
    )
    write(POOL, pool)


def consolidate_other_duplicates() -> None:
    # The strict per-cell-type projection fully replaces the older permissive path.
    old_projection = ROOT / "R/step_layer1_common_dictionary.R"
    if not old_projection.exists():
        raise RuntimeError("Expected obsolete Layer 1 projection file")
    old_projection.unlink()
    description = read("DESCRIPTION")
    line = "    'step_layer1_common_dictionary.R'\n"
    if description.count(line) != 1:
        raise RuntimeError("Layer 1 projection Collate entry not found exactly once")
    write("DESCRIPTION", description.replace(line, ""))

    # Keep one exported plotting API and call the boxplot implementation explicitly.
    rename_function(
        "R/condition_plot.R",
        "rc_plot_condition_reaction",
        ".rc_plot_condition_reaction_boxplot_core",
    )
    condition_plot = read("R/condition_plot.R")
    condition_plot = condition_plot.replace("#' @export\n", "#' @name rc_plot_condition_reaction\n", 1)
    write("R/condition_plot.R", condition_plot)
    violin = read("R/condition_plot_violin.R")
    alias = ".rc_plot_condition_reaction_boxplot_impl <- rc_plot_condition_reaction\n\n"
    if violin.count(alias) != 1:
        raise RuntimeError("Condition plot implementation alias not found")
    violin = violin.replace(alias, "")
    violin = violin.replace(
        ".rc_plot_condition_reaction_boxplot_impl(",
        ".rc_plot_condition_reaction_boxplot_core(",
    )
    write("R/condition_plot_violin.R", violin)

    remove_function("R/gpr_capacity.R", "rc_safe_scale")


def update_contract_tests() -> None:
    replacements = {
        "tests/testthat/test-layer2-corda-original-code-contract.R": {
            "body(RegCompassR:::.rc_corda_parent)":
                "body(RegCompassR:::.rc_corda_parent_core)",
            "body(RegCompassR:::.rc_complete_celltype_medium_corda_gem)":
                "body(RegCompassR:::.rc_complete_celltype_medium_corda_gem_core)",
        },
        "tests/testthat/test-layer2-corda-like.R": {
            "body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache)":
                "body(RegCompassR:::.rc_build_celltype_medium_union_gem_cache_core)",
        },
    }
    for path, values in replacements.items():
        text = read(path)
        for old, new in values.items():
            if old not in text:
                raise RuntimeError(f"Expected contract expression {old} in {path}")
            text = text.replace(old, new)
        write(path, text)


def main() -> None:
    consolidate_pool_wrappers()
    consolidate_other_duplicates()
    update_contract_tests()


if __name__ == "__main__":
    main()
