from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "R"


def require(text: str, marker: str, label: str) -> int:
    pos = text.find(marker)
    if pos < 0:
        raise RuntimeError(f"Missing {label}: {marker}")
    return pos


zz_path = R / "zz_reaction_annotations.R"
zz = zz_path.read_text(encoding="utf-8")

nonempty = require(zz, ".rc_ra_nonempty <- function", "base helper start")
group_evidence = require(zz, ".rc_ra_group_evidence <- function", "old evidence start")
annotation_from_object = require(
    zz, ".rc_ra_annotation_from_object <- function", "annotation join start"
)
pairwise_evidence = require(
    zz, ".rc_ra_pairwise_evidence <- function", "old pairwise evidence start"
)
group_evidence_join = require(
    zz, ".rc_ra_group_evidence_join <- function", "group evidence join start"
)
omnibus_evidence = require(
    zz, ".rc_ra_omnibus_evidence <- function", "old omnibus evidence start"
)
enrich_statistics = require(
    zz, ".rc_ra_enrich_statistics <- function", "statistics enrichment start"
)
stage_wrapper = require(zz, "# Enhanced Stage 6 assembly", "Stage 6 wrapper start")
select_start = require(
    zz, "#' Select scored reactions by metabolic genes", "gene selection start"
)
plot_caption = require(zz, ".rc_ra_plot_caption <- function", "plot caption start")
plot_wrapper = require(
    zz, "# Add biological labels to the existing single-reaction plot",
    "single plot wrapper start",
)
precomputed_plot = require(
    zz, ".rc_ra_plot_one_precomputed <- function", "batch plot helper start"
)

annotation_helper = r'''

.rc_ra_annotate_condition_plot <- function(
    plot, statistics, reaction_id, cell_type,
    target_direction = NULL, medium_scenario = NULL, title = NULL) {
  annotation_row <- data.frame()
  if (!is.list(statistics) || !is.data.frame(statistics$pairwise)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  required <- c(
    "reaction_id", "cell_type", "target_direction", "medium_scenario",
    "reaction_name", "tested_formula", "genes", "evidence_comparison"
  )
  if (!all(required %in% colnames(statistics$pairwise))) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }

  keep <- as.character(statistics$pairwise$reaction_id) == reaction_id &
    as.character(statistics$pairwise$cell_type) == cell_type
  if (!is.null(target_direction)) {
    keep <- keep &
      as.character(statistics$pairwise$target_direction) == target_direction
  }
  if (!is.null(medium_scenario)) {
    keep <- keep &
      as.character(statistics$pairwise$medium_scenario) == medium_scenario
  }
  annotation_row <- statistics$pairwise[keep, , drop = FALSE]
  if (!nrow(annotation_row)) {
    attr(plot, "reaction_annotation") <- annotation_row
    return(plot)
  }
  annotation_row <- annotation_row[1L, , drop = FALSE]
  reaction_name <- as.character(annotation_row$reaction_name[[1L]])
  evidence_text <- as.character(annotation_row$evidence_comparison[[1L]])
  if (is.null(title) && length(reaction_name) == 1L &&
      .rc_ra_nonempty(reaction_name)) {
    plot <- plot + ggplot2::labs(
      title = paste0(
        reaction_name, " (", reaction_id, ", ",
        annotation_row$target_direction[[1L]], ") in ", cell_type
      )
    )
  }
  caption <- .rc_ra_plot_caption(annotation_row, evidence_text)
  if (!is.null(caption) && length(caption) == 1L && !is.na(caption)) {
    plot <- plot + ggplot2::labs(caption = caption) + ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
  }
  attr(plot, "reaction_annotation") <- annotation_row
  plot
}
'''

reaction_annotations = (
    "# Biological reaction annotations and evidence joins.\n\n"
    + zz[nonempty:group_evidence]
    + zz[annotation_from_object:pairwise_evidence]
    + zz[group_evidence_join:omnibus_evidence]
    + zz[enrich_statistics:stage_wrapper]
    + zz[select_start:plot_caption]
    + zz[plot_caption:plot_wrapper]
    + annotation_helper
)
(R / "reaction_annotations.R").write_text(
    reaction_annotations, encoding="utf-8"
)

reaction_gene_plots = (
    "# Gene-centered collections of condition-reaction plots.\n\n"
    + zz[precomputed_plot:]
)
(R / "reaction_gene_plots.R").write_text(
    reaction_gene_plots, encoding="utf-8"
)

# Integrate annotation enrichment directly into condition statistics.
condition_stats_path = R / "condition_statistics.R"
condition_stats = condition_stats_path.read_text(encoding="utf-8")
marker = '  class(answer) <- c("regcompass_condition_statistics", "list")\n\n  if (!is.null(outdir)) {'
replacement = '''  class(answer) <- c("regcompass_condition_statistics", "list")
  annotation <- .rc_ra_annotation_from_object(x)
  answer <- .rc_ra_enrich_statistics(answer, annotation)

  if (!is.null(outdir)) {'''
if marker not in condition_stats:
    raise RuntimeError("Could not locate condition-statistics enrichment point")
condition_stats = condition_stats.replace(marker, replacement, 1)

save_marker = '''    saveRDS(answer, file.path(outdir, "condition_reaction_statistics.rds"))'''
save_replacement = '''    if (!is.null(annotation)) {
      .rc_write_tsv_gz(
        answer$reaction_catalog,
        file.path(outdir, "condition_reaction_catalog.tsv.gz")
      )
      .rc_write_tsv_gz(
        answer$reaction_evidence,
        file.path(outdir, "condition_reaction_evidence.tsv.gz")
      )
    }
    saveRDS(answer, file.path(outdir, "condition_reaction_statistics.rds"))'''
if save_marker not in condition_stats:
    raise RuntimeError("Could not locate condition-statistics save point")
condition_stats = condition_stats.replace(save_marker, save_replacement, 1)
condition_stats_path.write_text(condition_stats, encoding="utf-8")

# Integrate biological labels directly into the canonical single-reaction plot.
condition_plot_path = R / "condition_plot.R"
condition_plot = condition_plot_path.read_text(encoding="utf-8")
plot_marker = '''  attr(plot, "condition_statistics") <- statistics
  attr(plot, "plot_data") <- plot_data'''
plot_replacement = '''  plot <- .rc_ra_annotate_condition_plot(
    plot = plot,
    statistics = statistics,
    reaction_id = reaction_id,
    cell_type = cell_type,
    target_direction = target_direction,
    medium_scenario = medium_scenario,
    title = title
  )
  attr(plot, "condition_statistics") <- statistics
  attr(plot, "plot_data") <- plot_data'''
if plot_marker not in condition_plot:
    raise RuntimeError("Could not locate condition-plot annotation point")
condition_plot = condition_plot.replace(plot_marker, plot_replacement, 1)
condition_plot_path.write_text(condition_plot, encoding="utf-8")

# Integrate reaction annotation directly into Stage 6 assembly.
stepwise_path = R / "stepwise_workflow.R"
stepwise = stepwise_path.read_text(encoding="utf-8")
stage_marker = '''  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))'''
stage_replacement = '''  )
  result <- .rc_ra_attach_to_result(
    result = result,
    gem = gem,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_write_tsv_gz(
    result$reaction_catalog,
    file.path(outdir, "reaction_catalog.tsv.gz")
  )
  .rc_write_tsv_gz(
    result$reaction_evidence,
    file.path(outdir, "reaction_evidence_by_condition_celltype.tsv.gz")
  )
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))'''
if stage_marker not in stepwise:
    raise RuntimeError("Could not locate Stage 6 annotation point")
stepwise = stepwise.replace(stage_marker, stage_replacement, 1)
stepwise_path.write_text(stepwise, encoding="utf-8")

# Replace Collate entries with direct modules in dependency order.
description_path = ROOT / "DESCRIPTION"
description = description_path.read_text(encoding="utf-8")
old_collate = '''    'condition_layer1.R'
    'condition_statistics.R'
    'condition_plot.R'
    'layer1_builder.R'
    'parallel.R'
    'v170_microcompass_contract.R'
    'penalty.R'
    'reaction_roles.R'
    'stepwise_workflow.R'
    'zz_reaction_annotations.R'
    'reaction_evidence.R'
    'reaction_plot_annotations.R'
    'regcompass.R'
'''
new_collate = '''    'condition_layer1.R'
    'reaction_annotations.R'
    'reaction_evidence.R'
    'condition_statistics.R'
    'condition_plot.R'
    'reaction_gene_plots.R'
    'layer1_builder.R'
    'parallel.R'
    'v170_microcompass_contract.R'
    'penalty.R'
    'reaction_roles.R'
    'stepwise_workflow.R'
    'regcompass.R'
'''
if old_collate not in description:
    raise RuntimeError("Could not locate old reaction Collate block")
description = description.replace(old_collate, new_collate, 1)
description_path.write_text(description, encoding="utf-8")

# Update source-architecture tests to assert direct integration.
public_api_path = ROOT / "tests" / "testthat" / "test_public_api.R"
public_api = public_api_path.read_text(encoding="utf-8")
public_api = public_api.replace(
    '  expect_match(collate, "zz_reaction_annotations.R", fixed = TRUE)\n'
    '  expect_match(collate, "reaction_evidence.R", fixed = TRUE)',
    '  expect_match(collate, "reaction_annotations.R", fixed = TRUE)\n'
    '  expect_match(collate, "reaction_evidence.R", fixed = TRUE)\n'
    '  expect_match(collate, "reaction_gene_plots.R", fixed = TRUE)'
)
public_api_path.write_text(public_api, encoding="utf-8")

stage_contract_path = ROOT / "tests" / "testthat" / "test-stage-io-contracts.R"
stage_contract = stage_contract_path.read_text(encoding="utf-8")
old_stage_test = '''test_that("final results retain modules and add reaction interpretation", {
  assembly_text <- paste(deparse(body(.rc_step_results_core)), collapse = "\\n")
  expect_match(assembly_text, "condition_grn_meta_modules", fixed = TRUE)
  expect_match(assembly_text, "global_grn_meta_modules", fixed = TRUE)
  expect_match(assembly_text, "grn_metacell_group_coverage", fixed = TRUE)

  annotation_text <- paste(
    deparse(body(rc_regcompass_step_results)), collapse = "\\n"
  )
  expect_match(annotation_text, "reaction_catalog", fixed = TRUE)
  expect_match(annotation_text, "reaction_evidence", fixed = TRUE)
})'''
new_stage_test = '''test_that("final results retain modules and add reaction interpretation", {
  body_text <- paste(deparse(body(rc_regcompass_step_results)), collapse = "\\n")
  expect_match(body_text, "condition_grn_meta_modules", fixed = TRUE)
  expect_match(body_text, "global_grn_meta_modules", fixed = TRUE)
  expect_match(body_text, "grn_metacell_group_coverage", fixed = TRUE)
  expect_match(body_text, "reaction_catalog", fixed = TRUE)
  expect_match(body_text, "reaction_evidence", fixed = TRUE)
})'''
if old_stage_test not in stage_contract:
    raise RuntimeError("Could not locate wrapper-based Stage 6 test")
stage_contract = stage_contract.replace(old_stage_test, new_stage_test, 1)
stage_contract_path.write_text(stage_contract, encoding="utf-8")

rd_path = ROOT / "man" / "reaction_annotations.Rd"
rd = rd_path.read_text(encoding="utf-8")
rd = rd.replace(
    "% Please edit documentation in R/zz_reaction_annotations.R and R/zzz_reaction_evidence_hardening.R",
    "% Please edit documentation in R/reaction_annotations.R and R/reaction_evidence.R"
)
rd_path.write_text(rd, encoding="utf-8")

# Remove load-order wrapper sources after their reusable components are moved.
for obsolete in (
    R / "zz_reaction_annotations.R",
    R / "reaction_plot_annotations.R",
):
    obsolete.unlink()

print("Reaction annotations integrated directly into canonical source files.")
