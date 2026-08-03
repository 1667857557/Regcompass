from pathlib import Path

engine_path = Path('R/microcompass_engine.R')
engine = engine_path.read_text()

if '.rc_run_microcompass_engine <- function(' not in engine:
    marker = "#' Run directional minimum-evidence-discordance LPs\n"
    if marker not in engine:
        raise RuntimeError('microCOMPASS public-engine insertion marker missing')
    dispatcher = r'''.rc_run_microcompass_engine <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    mode = c("full_gem", "meta_module_gem"),
    reaction_membership = NULL, core_reactions = NULL,
    unit = c("metacell", "sample_celltype"),
    condition_col = "condition", sample_col = NULL,
    celltype_col = "cell_type", model_params = list(),
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    parallel = TRUE,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    BPPARAM = NULL,
    model_cache_override = NULL) {
  mode <- match.arg(mode)
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  if (identical(mode, "meta_module_gem")) {
    return(.rc_run_celltype_microcompass_engine(
      layer1 = layer1,
      gem = gem,
      target_reactions = target_reactions,
      medium_table = medium_table,
      medium_scenarios = medium_scenarios,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      unit = unit,
      condition_col = condition_col,
      sample_col = sample_col,
      celltype_col = celltype_col,
      model_params = model_params,
      omega = omega,
      target_direction = target_direction,
      parallel = parallel,
      solver = solver,
      flux_threshold = flux_threshold,
      BPPARAM = BPPARAM,
      model_cache_override = model_cache_override
    ))
  }
  .rc_run_shared_full_gem_engine(
    layer1 = layer1,
    gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = "full_gem",
    reaction_membership = NULL,
    core_reactions = NULL,
    unit = unit,
    condition_col = condition_col,
    sample_col = sample_col,
    celltype_col = celltype_col,
    model_params = model_params,
    omega = omega,
    target_direction = target_direction,
    parallel = parallel,
    solver = solver,
    flux_threshold = flux_threshold,
    BPPARAM = BPPARAM,
    model_cache_override = model_cache_override
  )
}

'''
    engine = engine.replace(marker, dispatcher + marker, 1)

engine = engine.replace(
    '          module_id = NA_character_,\nreaction_id = entry$reaction_id,\n',
    '          module_id = NA_character_,\n'
    '          reaction_id = entry$reaction_id,\n',
    1
)
engine = engine.replace(
    '      shared_gem_scope =\n'
    '        "one_full_gem_per_medium_shared_across_all_units",\n'
    'parallel_task = "shared_model_by_metacell_step2",\n',
    '      shared_gem_scope =\n'
    '        "one_full_gem_per_medium_shared_across_all_units",\n'
    '      parallel_task = "shared_model_by_metacell_step2",\n',
    1
)
engine_path.write_text(engine)

celltype_path = Path('R/celltype_microcompass_engine.R')
celltype = celltype_path.read_text()
old_coverage = '''  if (!setequal(cache_celltypes, unit_celltypes)) {
    stop(
      "Cell-type union GEMs and Layer 1 units cover different cell types.",
      call. = FALSE
    )
  }
'''
new_coverage = '''  if (is.null(model_cache_override)) {
    if (!setequal(cache_celltypes, unit_celltypes)) {
      stop(
        "Stage 5 union GEMs and Layer 1 units cover different cell types.",
        call. = FALSE
      )
    }
  } else if (!all(cache_celltypes %in% unit_celltypes)) {
    stop(
      "A reused union-GEM cache contains cell types absent from Layer 1.",
      call. = FALSE
    )
  }
'''
if celltype.count(old_coverage) != 1:
    raise RuntimeError('cell-type cache coverage contract was not found once')
celltype = celltype.replace(old_coverage, new_coverage, 1)
celltype_path.write_text(celltype)

test_path = Path('tests/testthat/test_union_gem_architecture.R')
test = test_path.read_text()
addition = r'''

test_that("microCOMPASS dispatches structural modes explicitly", {
  expect_true(exists(
    ".rc_run_microcompass_engine",
    envir = asNamespace("RegCompassR"),
    inherits = FALSE
  ))
  body_text <- paste(
    deparse(body(.rc_run_microcompass_engine)), collapse = "\n"
  )
  expect_true(grepl(
    ".rc_run_celltype_microcompass_engine", body_text, fixed = TRUE
  ))
  expect_true(grepl(
    ".rc_run_shared_full_gem_engine", body_text, fixed = TRUE
  ))
})

test_that("targeted cache reuse may cover a cell-type subset", {
  body_text <- paste(
    deparse(body(.rc_run_celltype_microcompass_engine)), collapse = "\n"
  )
  expect_true(grepl("is.null(model_cache_override)", body_text, fixed = TRUE))
  expect_true(grepl(
    "all(cache_celltypes %in% unit_celltypes)", body_text, fixed = TRUE
  ))
})
'''
if 'microCOMPASS dispatches structural modes explicitly' not in test:
    test_path.write_text(test.rstrip() + addition)

# The public execution path must exist exactly once and retired builders must
# remain absent.
all_r = '\n'.join(path.read_text() for path in Path('R').glob('*.R'))
if all_r.count('.rc_run_microcompass_engine <- function(') != 1:
    raise RuntimeError('microCOMPASS dispatcher must be defined exactly once')
for retired in (
    '.rc_build_medium_specific_union_gem_cache',
    '.rc_complete_medium_union_gem',
    '.rc_read_cached_union_gem'
):
    if retired in all_r:
        raise RuntimeError(f'retired global union implementation remains: {retired}')
