from pathlib import Path
import re

ROOT = Path('.')
R = ROOT / 'R'


def function_span(text, name):
    pattern = re.compile(rf'(?m)^{re.escape(name)}\s*<-\s*function\s*\(')
    match = pattern.search(text)
    if match is None:
        raise RuntimeError(f'Function not found: {name}')
    brace = text.find('{', match.end())
    if brace < 0:
        raise RuntimeError(f'Opening brace not found: {name}')
    depth = 0
    quote = None
    escaped = False
    comment = False
    i = brace
    while i < len(text):
        ch = text[i]
        if comment:
            if ch == '\n':
                comment = False
            i += 1
            continue
        if quote is not None:
            if escaped:
                escaped = False
            elif ch == '\\':
                escaped = True
            elif ch == quote:
                quote = None
            i += 1
            continue
        if ch in ('"', "'"):
            quote = ch
        elif ch == '#':
            comment = True
        elif ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = i + 1
                while end < len(text) and text[end] in ' \t':
                    end += 1
                while end < len(text) and text[end] == '\n':
                    end += 1
                return match.start(), end
        i += 1
    raise RuntimeError(f'Function did not close: {name}')


def replace_function(path, name, replacement):
    text = path.read_text()
    start, end = function_span(text, name)
    path.write_text(text[:start] + replacement.rstrip() + '\n\n' + text[end:])


def replace_once(path, old, new):
    text = path.read_text()
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f'{path}: expected one match, found {count}: {old[:80]}')
    path.write_text(text.replace(old, new, 1))


(R / 'meta_module_merge.R').write_text(r'''# Merge biological meta-modules within each cell type only.

.rc_merge_meta_modules_by_cell_type <- function(
    condition_modules,
    celltype_col = "cell_type",
    condition_col = "condition") {
  if (!is.list(condition_modules)) {
    stop("`condition_modules` must be a list.", call. = FALSE)
  }
  valid_name <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (!valid_name(celltype_col) || !valid_name(condition_col)) {
    stop("Cell-type and condition column names must be non-empty strings.",
         call. = FALSE)
  }

  names_to_merge <- c(
    "condition_fit_status", "tf_peak_gene_condition_all",
    "tf_peak_gene_condition", "supported_metabolic_genes",
    "core_gene_reaction", "reaction_membership", "meta_module_summary"
  )
  tables <- lapply(names_to_merge, function(name) {
    value <- condition_modules[[name]]
    if (is.data.frame(value)) value else data.frame()
  })
  names(tables) <- names_to_merge

  normalize_scope <- function(tab, label) {
    if (!is.data.frame(tab) || !nrow(tab) ||
        !all(c(celltype_col, "reaction_id") %in% colnames(tab))) {
      stop(
        "`", label, "` must contain non-empty `", celltype_col,
        "` and `reaction_id` columns.",
        call. = FALSE
      )
    }
    tab[[celltype_col]] <- trimws(as.character(tab[[celltype_col]]))
    tab$reaction_id <- trimws(as.character(tab$reaction_id))
    if (anyNA(tab[[celltype_col]]) || any(!nzchar(tab[[celltype_col]])) ||
        anyNA(tab$reaction_id) || any(!nzchar(tab$reaction_id))) {
      stop("Meta-module cell types and reaction IDs must be complete.",
           call. = FALSE)
    }
    tab
  }

  core <- normalize_scope(tables$core_gene_reaction, "core_gene_reaction")
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  biological <- normalize_scope(
    tables$reaction_membership, "reaction_membership"
  )
  cell_types <- sort(unique(as.character(biological[[celltype_col]])))
  core_cell_types <- unique(as.character(core[[celltype_col]]))
  missing_core <- setdiff(cell_types, core_cell_types)
  if (length(missing_core)) {
    stop(
      "No complete-GPR core reactions were produced for cell types: ",
      paste(missing_core, collapse = ", "),
      call. = FALSE
    )
  }

  subset_table <- function(tab, cell_type) {
    if (!is.data.frame(tab) || !nrow(tab)) return(tab)
    if (!celltype_col %in% colnames(tab)) return(tab[0, , drop = FALSE])
    value <- trimws(as.character(tab[[celltype_col]]))
    tab[value == cell_type, , drop = FALSE]
  }

  build_one <- function(cell_type) {
    one <- lapply(tables, subset_table, cell_type = cell_type)
    one_core <- core[core[[celltype_col]] == cell_type, , drop = FALSE]
    one_biological <- biological[
      biological[[celltype_col]] == cell_type, , drop = FALSE
    ]
    core_ids <- unique(as.character(one_core$reaction_id))
    biological_ids <- unique(as.character(one_biological$reaction_id))
    if (!length(core_ids) || !length(biological_ids)) {
      stop("Cell type `", cell_type, "` has an empty meta-module catalogue.",
           call. = FALSE)
    }
    missing <- setdiff(core_ids, biological_ids)
    if (length(missing)) {
      stop(
        "Cell type `", cell_type,
        "` has core reactions missing from its biological membership: ",
        paste(utils::head(missing, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    catalogue_id <- paste0(
      "CELLTYPE_META_MODULES::",
      utils::URLencode(cell_type, reserved = TRUE)
    )
    scoped_core <- data.frame(
      catalogue_id = catalogue_id,
      cell_type = cell_type,
      reaction_id = core_ids,
      is_core = TRUE,
      stringsAsFactors = FALSE
    )
    names(scoped_core)[names(scoped_core) == "cell_type"] <- celltype_col
    scoped_membership <- data.frame(
      catalogue_id = catalogue_id,
      cell_type = cell_type,
      reaction_id = biological_ids,
      is_core = biological_ids %in% core_ids,
      inclusion_stage = ifelse(
        biological_ids %in% core_ids,
        "celltype_union_core",
        "celltype_union_biological_member"
      ),
      stringsAsFactors = FALSE
    )
    names(scoped_membership)[names(scoped_membership) == "cell_type"] <-
      celltype_col

    status <- one$condition_fit_status
    source_groups <- if (
      is.data.frame(status) && "group_id" %in% colnames(status)
    ) {
      unique(as.character(status$group_id))
    } else {
      character()
    }
    source_conditions <- if (
      is.data.frame(status) && condition_col %in% colnames(status)
    ) {
      value <- trimws(as.character(status[[condition_col]]))
      sort(unique(value[!is.na(value) & nzchar(value)]))
    } else {
      character()
    }

    one$biological_reaction_membership <- one_biological
    one$merged_core_reactions <- scoped_core
    one$merged_reaction_membership <- scoped_membership
    one$schema_version <- "regcompass_celltype_meta_modules_v1"
    one$cell_type <- cell_type
    one$celltype_col <- celltype_col
    one$condition_col <- condition_col
    one$source_group_ids <- source_groups
    one$source_conditions <- source_conditions
    one$core_definition <-
      "condition_specific_complete_gpr_cores_unioned_within_cell_type"
    one$merge_source <-
      "condition_meta_module_reactions_deduplicated_within_cell_type"
    one$merge_scope <- "cell_type"
    one$is_gem <- FALSE
    one$fastcore_applied <- FALSE
    one
  }

  catalogues <- stats::setNames(lapply(cell_types, build_one), cell_types)
  bind <- function(field) {
    .rc_bind_frames_fill(lapply(catalogues, `[[`, field))
  }
  out <- tables
  out$biological_reaction_membership <- biological
  out$cell_type_catalogues <- catalogues
  out$merged_core_reactions <- bind("merged_core_reactions")
  out$merged_reaction_membership <- bind("merged_reaction_membership")
  out$schema_version <- "regcompass_celltype_merged_meta_modules_v1"
  out$celltype_col <- celltype_col
  out$condition_col <- condition_col
  out$cell_types <- cell_types
  out$core_definition <-
    "condition_specific_complete_gpr_cores_unioned_within_cell_type"
  out$merge_source <-
    "condition_meta_module_reactions_deduplicated_within_cell_type_only"
  out$merge_scope <- "cell_type"
  out$cross_celltype_merge <- FALSE
  out$is_gem <- FALSE
  out$fastcore_applied <- FALSE
  out
}
''')

celltype_cache = r'''# Build and validate one union GEM per cell type and medium.

.rc_safe_cache_token <- function(value) {
  paste(
    sprintf("%02x", as.integer(charToRaw(enc2utf8(as.character(value))))),
    collapse = ""
  )
}

.rc_validate_celltype_reaction_scope <- function(
    reaction_membership, core_reactions, celltype_col) {
  valid_name <- is.character(celltype_col) && length(celltype_col) == 1L &&
    !is.na(celltype_col) && nzchar(trimws(celltype_col))
  if (!valid_name) {
    stop("`celltype_col` must be one non-empty column name.", call. = FALSE)
  }
  validate <- function(tab, label) {
    required <- c(celltype_col, "reaction_id")
    if (!is.data.frame(tab) || !all(required %in% colnames(tab)) || !nrow(tab)) {
      stop("`", label, "` must contain cell type and reaction ID columns.",
           call. = FALSE)
    }
    tab[[celltype_col]] <- trimws(as.character(tab[[celltype_col]]))
    tab$reaction_id <- trimws(as.character(tab$reaction_id))
    if (anyNA(tab[[celltype_col]]) || any(!nzchar(tab[[celltype_col]])) ||
        anyNA(tab$reaction_id) || any(!nzchar(tab$reaction_id))) {
      stop("Cell-type-scoped reaction tables cannot contain missing values.",
           call. = FALSE)
    }
    tab
  }
  membership <- validate(reaction_membership, "reaction_membership")
  core <- validate(core_reactions, "core_reactions")
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  membership_types <- sort(unique(as.character(membership[[celltype_col]])))
  core_types <- sort(unique(as.character(core[[celltype_col]])))
  if (!setequal(membership_types, core_types)) {
    stop(
      "Reaction membership and core reactions must cover identical cell types.",
      call. = FALSE
    )
  }
  for (cell_type in membership_types) {
    member_ids <- unique(membership$reaction_id[
      membership[[celltype_col]] == cell_type
    ])
    core_ids <- unique(core$reaction_id[core[[celltype_col]] == cell_type])
    missing <- setdiff(core_ids, member_ids)
    if (length(missing)) {
      stop(
        "Cell type `", cell_type,
        "` has core reactions outside its own meta-module membership.",
        call. = FALSE
      )
    }
  }
  list(
    reaction_membership = membership,
    core_reactions = core,
    cell_types = membership_types
  )
}

.rc_celltype_target_reactions <- function(
    target_reactions, cell_type, celltype_col) {
  if (is.null(target_reactions)) return(NULL)
  if (is.data.frame(target_reactions)) {
    required <- c(celltype_col, "reaction_id")
    if (!all(required %in% colnames(target_reactions))) {
      stop("Scoped target reactions must contain cell type and reaction ID.",
           call. = FALSE)
    }
    value <- trimws(as.character(target_reactions[[celltype_col]]))
    return(unique(trimws(as.character(target_reactions$reaction_id[
      value == cell_type
    ]))))
  }
  unique(trimws(as.character(target_reactions)))
}

.rc_build_celltype_medium_union_gem_cache <- function(
    gem, reaction_membership, core_reactions,
    target_reactions = NULL, medium_scenarios = NULL,
    celltype_col = "cell_type",
    cache_dir = tempfile("RegCompassR_celltype_medium_union_gem_"),
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE) {
  target_direction <- match.arg(target_direction)
  scoped <- .rc_validate_celltype_reaction_scope(
    reaction_membership, core_reactions, celltype_col
  )
  reaction_membership <- scoped$reaction_membership
  core_reactions <- scoped$core_reactions
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("Cell-type union-GEM FASTCORE time limit must be positive.",
         call. = FALSE)
  }
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  scenarios <- scenarios[!is.na(scenarios) & nzchar(scenarios)]
  if (!length(scenarios)) {
    stop("No medium scenarios are available for union-GEM construction.",
         call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache <- list()
  summaries <- list()
  summary_index <- 0L
  for (cell_type in scoped$cell_types) {
    membership <- reaction_membership[
      reaction_membership[[celltype_col]] == cell_type, , drop = FALSE
    ]
    core <- core_reactions[
      core_reactions[[celltype_col]] == cell_type, , drop = FALSE
    ]
    selected_targets <- .rc_celltype_target_reactions(
      target_reactions, cell_type, celltype_col
    )
    if (!is.null(selected_targets)) {
      core <- core[core$reaction_id %in% selected_targets, , drop = FALSE]
    }
    if (!nrow(core)) {
      stop("No core targets remain for cell type `", cell_type, "`.",
           call. = FALSE)
    }

    for (scenario in scenarios) {
      medium <- medium_scenarios[
        as.character(medium_scenarios$medium_scenario_id) == scenario,
        , drop = FALSE
      ]
      if (!nrow(medium) ||
          (".no_constraints" %in% colnames(medium) &&
           all(medium$.no_constraints))) {
        medium <- NULL
      }
      model <- .rc_complete_celltype_medium_union_gem(
        gem = gem,
        reaction_membership = membership,
        core_reactions = core,
        cell_type = cell_type,
        medium_table = medium,
        target_direction = target_direction,
        solver = solver,
        time_limit = time_limit,
        fastcore_epsilon = fastcore_epsilon,
        max_support_reactions = max_support_reactions,
        strict = strict
      )
      model$shared_across_conditions <- TRUE
      model$shared_across_cell_types <- FALSE
      model$is_union_gem <- TRUE
      model$union_gem_cell_type <- cell_type
      model$union_gem_medium_scenario <- scenario
      model$union_gem_scope <-
        "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"

      file <- file.path(
        cache_dir,
        paste0(
          "union_gem__celltype_", .rc_safe_cache_token(cell_type),
          "__medium_", .rc_safe_cache_token(scenario), ".rds"
        )
      )
      saveRDS(model, file)
      checksum <- unname(tools::md5sum(file))
      summary_index <- summary_index + 1L
      summaries[[summary_index]] <- data.frame(
        cell_type = cell_type,
        medium_scenario = scenario,
        file = file,
        file_checksum = checksum,
        n_reactions = ncol(model$S),
        n_metabolites = nrow(model$S),
        n_celltype_biological_reactions =
          model$build_params$n_celltype_biological_reactions,
        n_celltype_fastcore_support_reactions =
          model$build_params$n_celltype_fastcore_support_reactions,
        target_status = model$target_status,
        build_strategy = "celltype_medium_union_gem",
        completion_stage =
          "celltype_specific_fastcore_after_condition_module_union",
        completion_time_limit = time_limit,
        stringsAsFactors = FALSE
      )

      if (!nrow(model$target_directions)) next
      for (i in seq_len(nrow(model$target_directions))) {
        reaction <- as.character(model$target_directions$reaction_id[[i]])
        direction <- as.character(model$target_directions$target_direction[[i]])
        key <- paste0(
          "celltype=", utils::URLencode(cell_type, reserved = TRUE),
          "::reaction=", utils::URLencode(reaction, reserved = TRUE),
          "::direction=", direction,
          "::medium=", utils::URLencode(scenario, reserved = TRUE)
        )
        cache[[key]] <- list(
          module_id = "CELLTYPE_MEDIUM_UNION_GEM",
          cell_type = cell_type,
          reaction_id = reaction,
          target_direction = direction,
          medium_scenario = scenario,
          condition = "all",
          file = file,
          file_checksum = checksum,
          build_strategy = "celltype_medium_union_gem"
        )
      }
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  attr(cache, "celltype_col") <- celltype_col
  attr(cache, "structural_scope") <- "cell_type_x_medium"
  cache
}

.rc_read_celltype_union_gem <- function(
    file, cell_type, medium_scenario, expected_checksum) {
  valid <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
  }
  if (!valid(file) || !file.exists(file)) {
    stop("A required cell-type union GEM cache file is unavailable.",
         call. = FALSE)
  }
  if (!valid(cell_type) || !valid(medium_scenario) ||
      !valid(expected_checksum)) {
    stop("Cell type, medium and checksum are required for cached union GEMs.",
         call. = FALSE)
  }
  observed_checksum <- unname(tools::md5sum(file))
  if (!identical(observed_checksum, as.character(expected_checksum))) {
    stop("A cell-type union GEM cache file failed its checksum check.",
         call. = FALSE)
  }
  model <- readRDS(file)
  if (!isTRUE(model$is_union_gem) ||
      !identical(as.character(model$union_gem_cell_type), cell_type) ||
      !identical(as.character(model$union_gem_medium_scenario), medium_scenario) ||
      !identical(
        as.character(model$union_gem_scope),
        "one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"
      )) {
    stop(
      "Cached union GEM provenance does not match its cell type and medium.",
      call. = FALSE
    )
  }
  model
}
'''
(R / 'celltype_union_gem_cache.R').write_text(celltype_cache)
old_cache = R / 'union_gem_cache.R'
if old_cache.exists():
    old_cache.unlink()

fastcore = R / 'fastcore.R'
text = fastcore.read_text()
text = text.replace(
    '.rc_complete_medium_union_gem <- function(\n    gem, reaction_membership, core_reactions, medium_table = NULL,',
    '.rc_complete_celltype_medium_union_gem <- function(\n    gem, reaction_membership, core_reactions, cell_type, medium_table = NULL,',
    1
)
text = text.replace(
    '  target_direction <- match.arg(target_direction)\n  if (!is.finite(fastcore_epsilon)',
    '  target_direction <- match.arg(target_direction)\n  if (!is.character(cell_type) || length(cell_type) != 1L ||\n      is.na(cell_type) || !nzchar(trimws(cell_type))) {\n    stop("`cell_type` must identify one non-empty cell type.", call. = FALSE)\n  }\n  cell_type <- trimws(cell_type)\n  if (!is.finite(fastcore_epsilon)',
    1
)
text = text.replace('"global_fastcore_completed"', '"celltype_fastcore_completed"')
text = text.replace(
    '"Global FASTCORE union-GEM completion failed for parent-feasible targets: ",',
    '"Cell-type FASTCORE union-GEM completion failed for parent-feasible targets in `",\n      cell_type, "`: ",'
)
text = text.replace('meta$global_fastcore_support <-', 'meta$celltype_fastcore_support <-')
text = text.replace(
    'meta$support_only <- meta$global_fastcore_support &',
    'meta$support_only <- meta$celltype_fastcore_support &'
)
text = text.replace('final$sample_id <- "global"', 'final$sample_id <- cell_type')
text = text.replace(
    'final$grn_module_id <- "MEDIUM_UNION_GEM"',
    'final$grn_module_id <- paste0("CELLTYPE_MEDIUM_UNION_GEM::", cell_type)\n  final$cell_type <- cell_type'
)
text = text.replace(
    '"one_medium_shared_across_conditions_and_metacells"',
    '"one_cell_type_one_medium_shared_across_conditions_and_matching_metacells"'
)
text = text.replace(
    'strategy = "medium_specific_union_gem",',
    'strategy = "celltype_medium_union_gem",\n    cell_type = cell_type,'
)
text = text.replace(
    'completion_stage = "single_global_fastcore_after_meta_module_merge",',
    'completion_stage =\n      "celltype_specific_fastcore_after_condition_module_union",'
)
text = text.replace(
    'n_merged_biological_reactions = length(biological),\n    n_global_fastcore_support_reactions = length(selected_support),',
    'n_celltype_biological_reactions = length(biological),\n    n_celltype_fastcore_support_reactions = length(selected_support),'
)
fastcore.write_text(text)

celltype_engine = r'''# Directional LP scoring on cell-type-specific structural models.

.rc_celltype_model_contract <- function(model_cache) {
  if (!is.list(model_cache) || !length(model_cache)) {
    stop("The cell-type model cache is empty.", call. = FALSE)
  }
  files <- vapply(model_cache, function(entry) as.character(entry$file), character(1))
  representative <- match(unique(files), files)
  rows <- lapply(representative, function(i) {
    entry <- model_cache[[i]]
    model <- .rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )
    validated <- rc_validate_gem(model)
    data.frame(
      cell_type = as.character(entry$cell_type),
      medium_scenario = as.character(entry$medium_scenario),
      model_file = as.character(entry$file),
      model_file_checksum = unname(tools::md5sum(entry$file)[[1L]]),
      n_reactions = ncol(validated$S),
      n_metabolites = nrow(validated$S),
      reaction_order_checksum =
        .rc_microcompass_object_checksum(colnames(validated$S)),
      metabolite_order_checksum =
        .rc_microcompass_object_checksum(rownames(validated$S)),
      stoichiometry_bounds_checksum =
        .rc_microcompass_object_checksum(list(
          S = validated$S, lb = validated$lb, ub = validated$ub
        )),
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      stringsAsFactors = FALSE
    )
  })
  answer <- do.call(rbind, rows)
  rownames(answer) <- NULL
  answer[order(answer$cell_type, answer$medium_scenario), , drop = FALSE]
}

.rc_run_celltype_microcompass_engine <- function(
    layer1, gem, target_reactions = NULL,
    medium_table = NULL, medium_scenarios = NULL,
    reaction_membership, core_reactions,
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
  unit <- match.arg(unit)
  solver <- match.arg(solver)
  target_direction <- match.arg(target_direction)
  medium_scenarios <- .rc_validate_shared_medium(
    medium_scenarios %||% medium_table
  )
  matrices <- rc_layer2_unit_matrices(
    layer1,
    if (identical(unit, "metacell")) "metacell" else "sample_celltype",
    sample_col, celltype_col, condition_col
  )
  unit_meta <- matrices$unit_meta
  if (!is.data.frame(unit_meta) || !celltype_col %in% colnames(unit_meta)) {
    stop("Layer 1 unit metadata lack the cell-type column.", call. = FALSE)
  }
  unit_id <- if ("unit_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$unit_id)
  } else if ("pool_id" %in% colnames(unit_meta)) {
    as.character(unit_meta$pool_id)
  } else {
    stop("Layer 1 unit metadata lack unit_id/pool_id.", call. = FALSE)
  }
  if (!setequal(unit_id, colnames(matrices$reaction_expression))) {
    stop("Layer 1 unit metadata and expression columns differ.", call. = FALSE)
  }
  unit_meta <- unit_meta[match(colnames(matrices$reaction_expression), unit_id),
                         , drop = FALSE]
  unit_celltype <- trimws(as.character(unit_meta[[celltype_col]]))
  if (anyNA(unit_celltype) || any(!nzchar(unit_celltype))) {
    stop("Layer 1 unit cell types must be complete.", call. = FALSE)
  }
  names(unit_celltype) <- colnames(matrices$reaction_expression)
  gem <- rc_annotate_reaction_roles(gem)

  if (!is.null(model_cache_override)) {
    if (!is.list(model_cache_override) || !length(model_cache_override) ||
        is.null(attr(model_cache_override, "summary"))) {
      stop("`model_cache_override` is not an audited cell-type model cache.",
           call. = FALSE)
    }
    model_cache <- model_cache_override
  } else {
    model_cache <- .rc_build_celltype_medium_union_gem_cache(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      celltype_col = celltype_col,
      cache_dir = model_params$cache_dir %||%
        tempfile("RegCompassR_celltype_union_gem_cache_"),
      target_direction = target_direction,
      solver = solver,
      time_limit = model_params$completion_time_limit %||% 300,
      fastcore_epsilon = model_params$fastcore_epsilon %||% 1e-4,
      max_support_reactions = model_params$max_support_reactions %||% 2000,
      strict = model_params$strict %||% TRUE
    )
  }
  if (!length(model_cache)) {
    stop("No cell-type union-GEM targets were available.", call. = FALSE)
  }
  cache_celltypes <- sort(unique(vapply(
    model_cache, function(entry) as.character(entry$cell_type), character(1)
  )))
  unit_celltypes <- sort(unique(unit_celltype))
  if (!setequal(cache_celltypes, unit_celltypes)) {
    stop(
      "Cell-type union GEMs and Layer 1 units cover different cell types.",
      call. = FALSE
    )
  }
  summary <- attr(model_cache, "summary")
  requested_media <- unique(as.character(medium_scenarios$medium_scenario_id))
  if (!is.data.frame(summary) ||
      !all(c("cell_type", "medium_scenario") %in% colnames(summary)) ||
      !setequal(unique(as.character(summary$medium_scenario)), requested_media)) {
    stop("Cell-type model cache and requested media differ.", call. = FALSE)
  }

  row_ids <- names(model_cache)
  units <- colnames(matrices$reaction_expression)
  model_keys <- vapply(model_cache, function(entry) entry$file, character(1))
  names(model_keys) <- row_ids
  unique_model_keys <- unique(unname(model_keys))
  representative_rows <- vapply(
    unique_model_keys,
    function(key) names(model_keys)[match(key, model_keys)],
    character(1)
  )
  names(representative_rows) <- unique_model_keys
  all_reactions <- unique(unlist(lapply(representative_rows, function(row_id) {
    entry <- model_cache[[row_id]]
    colnames(.rc_read_celltype_union_gem(
      entry$file, entry$cell_type, entry$medium_scenario,
      entry$file_checksum
    )$S)
  }), use.names = FALSE))
  penalties <- rc_compute_multiome_penalty(
    rc_align_reaction_expression(
      matrices$reaction_expression, all_reactions, NA_real_
    ),
    reaction_roles = gem$reaction_roles
  )
  vmax_cache <- .rc_build_microcompass_vmax_cache(
    model_cache = model_cache,
    mode = "meta_module_gem",
    model_keys = model_keys,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  vmax_cache_diagnostics <- do.call(rbind, lapply(row_ids, function(row_id) {
    entry <- model_cache[[row_id]]
    value <- vmax_cache[[row_id]]
    data.frame(
      row_id = row_id,
      cell_type = entry$cell_type,
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      vmax = as.numeric(value$vmax),
      feasible = isTRUE(value$feasible),
      status = as.character(value$status),
      computation_scope =
        "celltype_model_x_directional_target_once",
      stringsAsFactors = FALSE
    )
  }))
  rownames(vmax_cache_diagnostics) <- NULL

  penalty <- vmax <- matrix(
    NA_real_, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )
  feasible <- evaluated <- matrix(
    FALSE, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )
  task_rows <- lapply(unique_model_keys, function(model_key) {
    row_id <- representative_rows[[model_key]]
    cell_type <- as.character(model_cache[[row_id]]$cell_type)
    eligible <- names(unit_celltype)[unit_celltype == cell_type]
    if (!length(eligible)) {
      stop("No Layer 1 units match cell type `", cell_type, "`.",
           call. = FALSE)
    }
    data.frame(
      model_key = model_key,
      unit_id = eligible,
      stringsAsFactors = FALSE
    )
  })
  tasks <- do.call(rbind, task_rows)

  run_one_unit <- function(task) {
    model_key <- as.character(task$model_key)
    unit_id <- as.character(task$unit_id)
    selected_rows <- names(model_keys)[model_keys == model_key]
    first_entry <- model_cache[[selected_rows[[1L]]]]
    if (!identical(unit_celltype[[unit_id]], first_entry$cell_type)) {
      stop("A union GEM was assigned to a different cell type.", call. = FALSE)
    }
    model <- .rc_read_celltype_union_gem(
      first_entry$file, first_entry$cell_type,
      first_entry$medium_scenario, first_entry$file_checksum
    )
    target_results <- lapply(selected_rows, function(row_id) {
      entry <- model_cache[[row_id]]
      unit_penalty <- penalties$penalty[colnames(model$S), unit_id]
      target_index <- match(entry$reaction_id, colnames(model$S))
      if (is.na(target_index)) {
        stop("A target reaction is absent from its cell-type union GEM.",
             call. = FALSE)
      }
      target_penalty <- unit_penalty[[target_index]]
      evidence_available <- is.finite(unit_penalty)
      solver_penalty <- unit_penalty
      solver_penalty[!evidence_available] <- 0
      fit <- .rc_compass_step2_from_vmax_directional(
        S = model$S, lb = model$lb, ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = solver_penalty,
        vmax_result = vmax_cache[[row_id]],
        target_direction = entry$target_direction,
        omega = omega, solver = solver,
        flux_threshold = flux_threshold
      )
      target_available <- is.finite(target_penalty)
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = if (target_available) fit$penalty else NA_real_,
        vmax = fit$vmax,
        feasible = isTRUE(fit$feasible),
        evaluated = isTRUE(fit$feasible) && target_available,
        diagnostics = data.frame(
          row_id = row_id,
          unit_id = unit_id,
          module_id = "CELLTYPE_MEDIUM_UNION_GEM",
          cell_type = entry$cell_type,
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          medium_scenario = entry$medium_scenario,
          condition = "all",
          strict_feasible = isTRUE(fit$feasible),
          solver_status = fit$solver_status,
          step1_status = fit$step1_status,
          step2_status = fit$step2_status,
          target_status = model$target_status %||%
            if (isTRUE(fit$feasible)) "ok" else "structurally_infeasible",
          objective_value = if (target_available) fit$penalty else NA_real_,
          vmax = fit$vmax,
          vmax_reused_from_celltype_cache = TRUE,
          target_expression_available = target_available,
          objective_evidence_fraction = mean(evidence_available),
          unavailable_objective_terms = sum(!evidence_available),
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      results = target_results,
      diagnostics = do.call(rbind, lapply(target_results, `[[`, "diagnostics"))
    )
  }

  grouped <- rc_parallel_lapply(
    split(tasks, seq_len(nrow(tasks))),
    function(task) run_one_unit(task[1, , drop = FALSE]),
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  results <- unlist(lapply(grouped, `[[`, "results"), recursive = FALSE)
  for (result in results) {
    penalty[result$row_id, result$unit_id] <- result$penalty
    vmax[result$row_id, result$unit_id] <- result$vmax
    feasible[result$row_id, result$unit_id] <- result$feasible
    evaluated[result$row_id, result$unit_id] <- result$evaluated
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {
    data.frame(
      cell_type = entry$cell_type,
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      stringsAsFactors = FALSE
    )
  })))
  model_diagnostics <- .rc_bind_frames_fill(lapply(
    representative_rows, function(row_id) {
      entry <- model_cache[[row_id]]
      model <- .rc_read_celltype_union_gem(
        entry$file, entry$cell_type, entry$medium_scenario,
        entry$file_checksum
      )
      out <- model$closure_diagnostics %||% data.frame()
      if (nrow(out)) {
        out$cell_type <- entry$cell_type
        out$medium_scenario <- entry$medium_scenario
      }
      out
    }
  ))
  reuse <- table(unit_celltype)
  list(
    score = score,
    penalty = penalty,
    vmax = vmax,
    feasible = feasible,
    evaluated = evaluated,
    target_direction = directions,
    direction_diagnostics = directions,
    medium_scenarios = medium_scenarios,
    model_mode = "meta_module_gem",
    shared_model_cache = model_cache,
    model_cache_summary = summary,
    structural_model_contract = .rc_celltype_model_contract(model_cache),
    model_diagnostics = model_diagnostics,
    vmax_cache_diagnostics = vmax_cache_diagnostics,
    lp_diagnostics = .rc_bind_frames_fill(lapply(grouped, `[[`, "diagnostics")),
    penalty_components = penalties$components,
    evidence_policy = penalties$evidence_policy,
    evidence_policy_detail = penalties$evidence_policy_detail,
    unit_meta = unit_meta,
    params = list(
      unit = unit,
      omega = omega,
      target_direction = target_direction,
      shared_gem = TRUE,
      shared_across_conditions = TRUE,
      shared_across_cell_types = FALSE,
      shared_gem_scope =
        "one_union_gem_per_cell_type_per_medium_shared_within_cell_type",
      structural_scope = "cell_type_x_medium",
      parallel_task = "celltype_model_by_matching_metacell_step2",
      vmax_computation_scope =
        "celltype_model_x_directional_target_once",
      vmax_solve_count = length(vmax_cache),
      vmax_reuse_by_cell_type = stats::setNames(
        as.integer(reuse), names(reuse)
      ),
      flux_threshold = flux_threshold,
      scoring_time_limit = "none"
    ),
    method = paste(
      "microCOMPASS directional LP on cell-type-specific medium union GEMs",
      "after independent FASTCORE completion within each cell type"
    )
  )
}
'''
(R / 'celltype_microcompass_engine.R').write_text(celltype_engine)

# Add the cell-type field to row-ID parsing.
parser = r'''rc_parse_microcompass_row_id <- function(x) {
  x <- as.character(x)
  required <- c("reaction", "direction", "medium")

  parse_one <- function(id) {
    fields <- strsplit(id, "::", fixed = TRUE)[[1L]]
    equals <- regexpr("=", fields, fixed = TRUE)
    if (any(equals < 2L)) {
      stop(
        "microCOMPASS row IDs must use labeled key=value fields.",
        call. = FALSE
      )
    }
    keys <- substring(fields, 1L, equals - 1L)
    values <- utils::URLdecode(substring(fields, equals + 1L))
    if (anyDuplicated(keys)) {
      stop("microCOMPASS row-ID fields must be unique.", call. = FALSE)
    }
    required_counts <- table(factor(keys, levels = required))
    required_values <- values[match(required, keys)]
    invalid <- any(required_counts != 1L) ||
      anyNA(required_values) || any(!nzchar(trimws(required_values))) ||
      !required_values[[2L]] %in% c("forward", "reverse")
    if (invalid) {
      stop(
        paste(
          "microCOMPASS row IDs require one non-empty reaction, direction",
          "and medium field; direction must be forward or reverse."
        ),
        call. = FALSE
      )
    }
    named <- stats::setNames(values, keys)
    value <- function(name) {
      hit <- named[name]
      if (!length(hit) || is.na(hit[[1L]]) ||
          !nzchar(trimws(hit[[1L]]))) return(NA_character_)
      unname(hit[[1L]])
    }
    data.frame(
      sample_id = value("sample"),
      module_id = value("module"),
      cell_type = value("celltype"),
      reaction_id = value("reaction"),
      target_direction = value("direction"),
      medium_scenario = value("medium"),
      condition = value("condition"),
      stringsAsFactors = FALSE
    )
  }
  rows <- lapply(x, parse_one)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
'''
replace_function(R / 'microcompass.R', 'rc_parse_microcompass_row_id', parser)

# Load cell-type caches with exact scope validation.
loader = r'''.rc_load_microcompass_model <- function(entry, mode) {
  if (identical(mode, "meta_module_gem")) {
    return(.rc_read_celltype_union_gem(
      file = entry$file,
      cell_type = entry$cell_type,
      medium_scenario = entry$medium_scenario,
      expected_checksum = entry$file_checksum %||% NA_character_
    ))
  }
  .rc_cache_gem(entry)
}
'''
replace_function(R / 'microcompass_engine.R', '.rc_load_microcompass_model', loader)

# Rename the previous generic implementation and dispatch meta-module mode to
# the cell-type-specific engine. The old implementation remains the full-GEM
# execution engine and is never used to build union GEMs.
engine_path = R / 'microcompass_engine.R'
engine_text = engine_path.read_text()
engine_text = engine_text.replace(
    '.rc_run_microcompass_engine <- function(',
    '.rc_run_shared_full_gem_engine <- function(',
    1
)
marker = "#' Run directional minimum-evidence-discordance LPs\n"
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
      layer1 = layer1, gem = gem,
      target_reactions = target_reactions,
      medium_table = medium_table,
      medium_scenarios = medium_scenarios,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      unit = unit, condition_col = condition_col,
      sample_col = sample_col, celltype_col = celltype_col,
      model_params = model_params, omega = omega,
      target_direction = target_direction, parallel = parallel,
      solver = solver, flux_threshold = flux_threshold,
      BPPARAM = BPPARAM,
      model_cache_override = model_cache_override
    ))
  }
  .rc_run_shared_full_gem_engine(
    layer1 = layer1, gem = gem,
    target_reactions = target_reactions,
    medium_table = medium_table,
    medium_scenarios = medium_scenarios,
    mode = "full_gem",
    reaction_membership = NULL,
    core_reactions = NULL,
    unit = unit, condition_col = condition_col,
    sample_col = sample_col, celltype_col = celltype_col,
    model_params = model_params, omega = omega,
    target_direction = target_direction, parallel = parallel,
    solver = solver, flux_threshold = flux_threshold,
    BPPARAM = BPPARAM,
    model_cache_override = model_cache_override
  )
}

'''
if marker not in engine_text:
    raise RuntimeError('microcompass engine insertion marker not found')
engine_path.write_text(engine_text.replace(marker, dispatcher + marker, 1))

# Stage 3 calls the cell-type merge contract.
step3 = R / 'step_meta_modules.R'
replace_once(
    step3,
    '  merged_modules <- .rc_merge_meta_module_catalogue(condition_modules)\n',
    '  merged_modules <- .rc_merge_meta_modules_by_cell_type(\n'
    '    condition_modules,\n'
    '    celltype_col = metacells$params$celltype_col,\n'
    '    condition_col = metacells$params$condition_col\n'
    '  )\n'
)
replace_once(
    step3,
    '      merge_creates_gem = FALSE,\n',
    '      merge_creates_gem = FALSE,\n'
    '      meta_module_merge_scope = "cell_type",\n'
    '      cross_celltype_merge = FALSE,\n'
)

# Layer 2 uses scoped core rows and records the structural key.
step2 = R / 'step_layer2.R'
replace_once(
    step2,
    '  targets <- unique(as.character(\n    catalogue$merged_core_reactions$reaction_id\n  ))\n',
    '  required_catalogue_cols <- c(params$celltype_col, "reaction_id")\n'
    '  if (!is.list(catalogue$cell_type_catalogues) ||\n'
    '      !all(required_catalogue_cols %in%\n'
    '           colnames(catalogue$merged_core_reactions)) ||\n'
    '      !all(required_catalogue_cols %in%\n'
    '           colnames(catalogue$merged_reaction_membership)) ||\n'
    '      !identical(catalogue$merge_scope, "cell_type") ||\n'
    '      isTRUE(catalogue$cross_celltype_merge)) {\n'
    '    stop("Meta-modules are not partitioned by cell type.", call. = FALSE)\n'
    '  }\n'
    '  targets <- unique(as.character(\n'
    '    catalogue$merged_core_reactions$reaction_id\n'
    '  ))\n'
)
replace_once(
    step2,
    '    target_reactions = targets,\n',
    '    target_reactions = if (identical(model_mode, "meta_module_gem")) {\n'
    '      catalogue$merged_core_reactions\n'
    '    } else {\n'
    '      targets\n'
    '    },\n'
)
replace_once(
    step2,
    '    "one medium-specific union GEM; single global FASTCORE completion"\n',
    '    paste(\n'
    '      "one union GEM per cell type and medium; FASTCORE runs",\n'
    '      "independently within each cell type"\n'
    '    )\n'
)

# DESCRIPTION and Collate reflect the functional source names.
desc = ROOT / 'DESCRIPTION'
dtext = desc.read_text()
dtext = dtext.replace(
    'RNA support, regulatory\n    effects, GPR rules, shared union GEM construction and directional metabolic\n',
    'RNA support, regulatory\n    effects, GPR rules, cell-type-specific union GEM construction and directional\n'
)
dtext = dtext.replace("    'union_gem_cache.R'\n", "    'celltype_union_gem_cache.R'\n")
dtext = dtext.replace(
    "    'microcompass_engine.R'\n",
    "    'microcompass_engine.R'\n    'celltype_microcompass_engine.R'\n"
)
desc.write_text(dtext)

# Basic source-level invariants.
for obsolete in ('union_gem_cache.R',):
    if (R / obsolete).exists():
        raise RuntimeError(f'Obsolete source remains: {obsolete}')
all_text = '\n'.join(path.read_text() for path in R.glob('*.R'))
for retired in (
    '.rc_merge_meta_module_catalogue <- function',
    '.rc_build_medium_specific_union_gem_cache <- function',
    '.rc_complete_medium_union_gem <- function',
    '.rc_read_cached_union_gem <- function'
):
    if retired in all_text:
        raise RuntimeError(f'Retired global implementation remains: {retired}')
