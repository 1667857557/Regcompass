from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def write(path, text):
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected one match, found {count}")
    return text.replace(old, new, 1)


def regex_once(text, pattern, replacement, label):
    value, count = re.subn(pattern, replacement, text, count=1, flags=re.S)
    if count != 1:
        raise RuntimeError(f"{label}: expected one regex match, found {count}")
    return value


# Prefer compatible package-bundled GEMs, while retaining the original
# download/conversion path for updates and custom versions.
hg = read("R/humangem.R")
hg = replace_once(
    hg,
    "#' Downloads and converts a pinned official SysBioChalmers GEM release. Human\n#' mode uses Human-GEM and mouse mode uses Mouse-GEM directly; mouse genes are\n#' retained as mouse symbols and are not converted through human orthologues.\n",
    "#' Loads a bundled pinned SysBioChalmers GEM by default. The official download\n#' and conversion path remains available for custom or updated releases. Mouse\n#' genes remain mouse symbols and are not converted through human orthologues.\n",
    "prepare GEM description",
)
hg = replace_once(
    hg,
    "#' @param force_download Re-download and rebuild an existing cached model.\n#' @param allow_latest Permit the unpinned `version = \"latest\"` mode.\n",
    "#' @param force_download Re-download and rebuild an existing cached model.\n#' @param allow_latest Permit the unpinned `version = \"latest\"` mode.\n#' @param source Model source. `auto` uses a compatible cache, then the bundled\n#'   pinned model, then downloads only when required. `bundled` forbids network\n#'   fallback; `download` retains the original preparation path.\n",
    "prepare GEM source docs",
)
hg = replace_once(
    hg,
    "    save_rds = NULL,\n    force_download = FALSE,\n    allow_latest = FALSE) {\n  species <- match.arg(species)\n  spec <- .rc_species_gem_spec(species, version)\n",
    "    save_rds = NULL,\n    force_download = FALSE,\n    allow_latest = FALSE,\n    source = c(\"auto\", \"bundled\", \"download\")) {\n  species <- match.arg(species)\n  source <- match.arg(source)\n  spec <- .rc_species_gem_spec(species, version)\n",
    "prepare GEM signature",
)
hg = replace_once(
    hg,
    "  if (isTRUE(force_download) && file.exists(save_rds)) {\n    unlink(save_rds, force = TRUE)\n  }\n  cached <- if (isTRUE(force_download)) {\n    NULL\n  } else {\n    .rc_load_compatible_species_gem(save_rds, spec)\n  }\n  if (!is.null(cached)) return(cached)\n\n  ref <- if (identical(spec$version, \"latest\")) {\n",
    "  if (isTRUE(force_download) && identical(source, \"bundled\")) {\n    stop(\"`force_download = TRUE` cannot be combined with `source = 'bundled'`.\",\n         call. = FALSE)\n  }\n  if (isTRUE(force_download)) source <- \"download\"\n  if (isTRUE(force_download) && file.exists(save_rds)) {\n    unlink(save_rds, force = TRUE)\n  }\n  cached <- if (isTRUE(force_download) || identical(source, \"bundled\")) {\n    NULL\n  } else {\n    .rc_load_compatible_species_gem(save_rds, spec)\n  }\n  if (!is.null(cached)) return(cached)\n\n  if (!identical(source, \"download\")) {\n    bundled <- .rc_load_bundled_species_gem(spec)\n    if (!is.null(bundled)) return(bundled)\n    if (identical(source, \"bundled\")) {\n      stop(\n        \"No bundled \", spec$repository_name, \" model matches version `\",\n        spec$version, \"`.\", call. = FALSE\n      )\n    }\n  }\n\n  ref <- if (identical(spec$version, \"latest\")) {\n",
    "prepare GEM load order",
)
hg = replace_once(
    hg,
    "    force_download = FALSE,\n    allow_latest = FALSE) {\n  rc_prepare_gem(\n    species = \"human\",\n",
    "    force_download = FALSE,\n    allow_latest = FALSE,\n    source = c(\"auto\", \"bundled\", \"download\")) {\n  source <- match.arg(source)\n  rc_prepare_gem(\n    species = \"human\",\n",
    "human helper signature",
)
hg = replace_once(
    hg,
    "    force_download = force_download,\n    allow_latest = allow_latest\n  )\n}\n\n#' Prepare Mouse-GEM for RegCompass\n",
    "    force_download = force_download,\n    allow_latest = allow_latest,\n    source = source\n  )\n}\n\n#' Prepare Mouse-GEM for RegCompass\n",
    "human helper source",
)
hg = replace_once(
    hg,
    "    force_download = FALSE,\n    allow_latest = FALSE) {\n  rc_prepare_gem(\n    species = \"mouse\",\n",
    "    force_download = FALSE,\n    allow_latest = FALSE,\n    source = c(\"auto\", \"bundled\", \"download\")) {\n  source <- match.arg(source)\n  rc_prepare_gem(\n    species = \"mouse\",\n",
    "mouse helper signature",
)
hg = replace_once(
    hg,
    "    force_download = force_download,\n    allow_latest = allow_latest\n  )\n}\n\nrc_validate_species_gem <- function",
    "    force_download = force_download,\n    allow_latest = allow_latest,\n    source = source\n  )\n}\n\nrc_validate_species_gem <- function",
    "mouse helper source",
)
hg = replace_once(
    hg,
    "rc_download_species_gem <- function(\n",
    "#' Download and parse an official species GEM release\n#'\n#' Retained for users rebuilding the bundled files or preparing a newer pinned\n#' upstream release.\n#'\n#' @return Parsed official model tables with source archive metadata.\n#' @export\nrc_download_species_gem <- function(\n",
    "download export",
)
write("R/humangem.R", hg)


# Add a monitor to each public step. The computational bodies are unchanged;
# each step now returns one timing row and writes step_timing.tsv.
sw = read("R/stepwise_workflow.R")
signature_changes = [
    (
        "    pando_args = list(),\n    parallel = TRUE,\n    BPPARAM = NULL) {",
        "    pando_args = list(),\n    parallel = TRUE,\n    BPPARAM = NULL,\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "GRN progress argument",
    ),
    (
        "    fragment_files = FALSE,\n    metacell_args = list()) {",
        "    fragment_files = FALSE,\n    metacell_args = list(),\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "metacell progress argument",
    ),
    (
        "    grn, metacells, gem, outdir,\n    layer1_args = list()) {",
        "    grn, metacells, gem, outdir,\n    layer1_args = list(),\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "meta-module progress argument",
    ),
    (
        "    gene_half_saturation = getOption(\"RegCompassR.cpm_half_saturation\", 1),\n    parallel = TRUE,\n    BPPARAM = NULL) {",
        "    gene_half_saturation = getOption(\"RegCompassR.cpm_half_saturation\", 1),\n    parallel = TRUE,\n    BPPARAM = NULL,\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "Layer 1 progress argument",
    ),
    (
        "    model_mode = c(\"meta_module_gem\", \"full_gem\"),\n    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {",
        "    model_mode = c(\"meta_module_gem\", \"full_gem\"),\n    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "Layer 2 progress argument",
    ),
    (
        "    grn, metacells, meta_modules, layer1, layer2, gem, outdir,\n    species = c(\"auto\", \"human\", \"mouse\")) {",
        "    grn, metacells, meta_modules, layer1, layer2, gem, outdir,\n    species = c(\"auto\", \"human\", \"mouse\"),\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {",
        "results progress argument",
    ),
]
for old, new, label in signature_changes:
    sw = replace_once(sw, old, new, label)

for function_name, stage in [
    ("rc_regcompass_step_grn", "grn"),
    ("rc_regcompass_step_metacells", "metacells"),
    ("rc_regcompass_step_meta_modules", "meta_modules"),
    ("rc_regcompass_step_layer1", "layer1"),
    ("rc_regcompass_step_layer2", "layer2"),
    ("rc_regcompass_step_results", "results"),
]:
    pattern = rf"({function_name} <- function\([\s\S]*?\) \{{\n)"
    replacement = (
        rf"\1  monitor <- .rc_step_monitor_start(\"{stage}\", outdir, progress)\n"
        "  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)\n"
    )
    sw = regex_once(sw, pattern, replacement, f"{function_name} monitor")

finish_changes = [
    (
        "  class(answer) <- c(\"regcompass_grn_step\", \"list\")\n  saveRDS(answer, file.path(outdir, \"step_grn.rds\"))",
        "  class(answer) <- c(\"regcompass_grn_step\", \"list\")\n  answer <- .rc_step_monitor_finish(answer, monitor)\n  saveRDS(answer, file.path(outdir, \"step_grn.rds\"))",
        "finish GRN",
    ),
    (
        "  class(answer) <- c(\"regcompass_metacell_step\", \"list\")\n  saveRDS(answer, file.path(outdir, \"step_metacells.rds\"))",
        "  class(answer) <- c(\"regcompass_metacell_step\", \"list\")\n  answer <- .rc_step_monitor_finish(answer, monitor)\n  saveRDS(answer, file.path(outdir, \"step_metacells.rds\"))",
        "finish metacells",
    ),
    (
        "  class(answer) <- c(\"regcompass_meta_module_step\", \"list\")\n  saveRDS(condition_modules, file.path(outdir, \"condition_meta_modules.rds\"))",
        "  class(answer) <- c(\"regcompass_meta_module_step\", \"list\")\n  answer <- .rc_step_monitor_finish(answer, monitor)\n  saveRDS(condition_modules, file.path(outdir, \"condition_meta_modules.rds\"))",
        "finish meta-modules",
    ),
    (
        "  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)\n  saveRDS(layer1, file.path(outdir, \"step_layer1.rds\"))\n  layer1",
        "  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)\n  layer1 <- .rc_step_monitor_finish(layer1, monitor)\n  saveRDS(layer1, file.path(outdir, \"step_layer1.rds\"))\n  layer1",
        "finish Layer 1",
    ),
    (
        "  rc_export_microcompass(answer, outdir)\n  saveRDS(answer, file.path(outdir, \"step_layer2.rds\"))",
        "  answer <- .rc_step_monitor_finish(answer, monitor)\n  rc_export_microcompass(answer, outdir)\n  saveRDS(answer, file.path(outdir, \"step_layer2.rds\"))",
        "finish Layer 2",
    ),
    (
        "  saveRDS(comparison, file.path(outdir, \"step_comparison.rds\"))\n  saveRDS(result, file.path(outdir, \"regcompass_result.rds\"))",
        "  result <- .rc_step_monitor_finish(result, monitor)\n  saveRDS(comparison, file.path(outdir, \"step_comparison.rds\"))\n  saveRDS(result, file.path(outdir, \"regcompass_result.rds\"))",
        "finish results",
    ),
]
for old, new, label in finish_changes:
    sw = replace_once(sw, old, new, label)
sw = sw.replace('version = "1.8.2"', 'version = "1.8.3"')
write("R/stepwise_workflow.R", sw)


# Target-union is also a public analysis stage and receives the same monitor.
tu = read("R/target_union.R")
tu = replace_once(
    tu,
    "    gene_match = c(\"complete_gpr\", \"any_direct\"),\n    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {\n",
    "    gene_match = c(\"complete_gpr\", \"any_direct\"),\n    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,\n    progress = getOption(\"RegCompassR.progress\", TRUE)) {\n  monitor <- .rc_step_monitor_start(\"target_union\", outdir, progress)\n  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)\n",
    "target-union monitor",
)
tu = replace_once(
    tu,
    "  class(answer) <- c(\"regcompass_target_union_step\", \"list\")\n  saveRDS(answer, file.path(outdir, \"step_target_union.rds\"))",
    "  class(answer) <- c(\"regcompass_target_union_step\", \"list\")\n  answer <- .rc_step_monitor_finish(answer, monitor)\n  saveRDS(answer, file.path(outdir, \"step_target_union.rds\"))",
    "target-union finish",
)
write("R/target_union.R", tu)


# Keep all user-facing version assertions synchronized.
for path in [
    "README.md",
    "docs/functions.md",
    "docs/tutorial-01-quick-start.md",
    "docs/tutorial-02-stepwise-audit.md",
    "docs/tutorial-03-advanced-restart.md",
    "docs/workflow.md",
    "docs/run-modes-and-stepwise-workflow.md",
    "docs/stage-interface-contracts.md",
    "vignettes/regcompass-workflow.Rmd",
    "tests/testthat/test-stage-io-contracts.R",
    "tests/testthat/test-vignette-contract.R",
]:
    target = ROOT / path
    if target.exists():
        write(path, read(path).replace("1.8.2", "1.8.3"))

print("Applied RegCompassR 1.8.3 compatibility patches.")
