# Build and validate one union GEM per cell type and medium.

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
  context <- .rc_layer2_completion_context
  if (isTRUE(context$active) &&
      identical(context$model_completion, "corda2") &&
      .rc_is_corda2_options(context$corda_options)) {
    return(.rc_build_celltype_medium_corda_cache(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      target_reactions = target_reactions,
      medium_scenarios = medium_scenarios,
      celltype_col = celltype_col,
      cache_dir = cache_dir,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    ))
  }

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
