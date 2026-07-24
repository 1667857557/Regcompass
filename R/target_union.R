# Re-score annotation-related reactions in a previously built global union GEM.

.rc_target_union_normalize_ids <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

.rc_target_union_core_rows <- function(
    gem, available_core_reactions,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct")) {
  gene_match <- match.arg(gene_match)
  validated <- rc_validate_gem(gem)
  available <- intersect(
    .rc_target_union_normalize_ids(available_core_reactions),
    validated$reactions
  )
  if (!length(available)) {
    stop("The previous analysis contains no valid global core reactions.",
         call. = FALSE)
  }

  requested_reactions <- .rc_target_union_normalize_ids(core_reaction_ids)
  requested_genes <- toupper(.rc_target_union_normalize_ids(core_genes))
  if (!length(requested_reactions) && !length(requested_genes)) {
    stop(
      "Supply at least one `core_reaction_ids` or `core_genes` value.",
      call. = FALSE
    )
  }

  missing_reactions <- setdiff(requested_reactions, validated$reactions)
  if (length(missing_reactions)) {
    stop(
      "Selected reactions are absent from the GEM: ",
      paste(utils::head(missing_reactions, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  not_previous_core <- setdiff(requested_reactions, available)
  if (length(not_previous_core)) {
    stop(
      "Selected reaction IDs were not core reactions in the previous LP analysis: ",
      paste(utils::head(not_previous_core, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  gpr <- gem$gpr_table
  gene_reactions <- character()
  gene_source <- list()
  if (length(requested_genes)) {
    required <- c("reaction_id", "and_group_id", "gene")
    if (!is.data.frame(gpr) || !all(required %in% colnames(gpr))) {
      stop(
        paste(
          "Gene-selected cores require a GEM `gpr_table` containing",
          "reaction_id, and_group_id and gene."
        ),
        call. = FALSE
      )
    }
    gpr <- unique(gpr[, required, drop = FALSE])
    gpr$reaction_id <- trimws(as.character(gpr$reaction_id))
    gpr$and_group_id <- as.character(gpr$and_group_id)
    gpr$gene <- toupper(trimws(as.character(gpr$gene)))
    gpr <- gpr[
      gpr$reaction_id %in% validated$reactions &
        !is.na(gpr$gene) & nzchar(gpr$gene),
      , drop = FALSE
    ]
    absent_genes <- setdiff(requested_genes, unique(gpr$gene))
    if (length(absent_genes)) {
      stop(
        "Selected core genes do not map to GEM GPR rules: ",
        paste(utils::head(absent_genes, 10L), collapse = ", "),
        call. = FALSE
      )
    }

    if (identical(gene_match, "any_direct")) {
      mapped <- unique(gpr$reaction_id[gpr$gene %in% requested_genes])
    } else {
      group_key <- paste(gpr$reaction_id, gpr$and_group_id, sep = "\001")
      groups <- split(seq_len(nrow(gpr)), group_key)
      complete <- vapply(groups, function(rows) {
        all(unique(gpr$gene[rows]) %in% requested_genes)
      }, logical(1))
      mapped <- unique(vapply(groups[complete], function(rows) {
        gpr$reaction_id[rows[[1L]]]
      }, character(1)))
    }
    gene_reactions <- intersect(mapped, available)
    if (!length(gene_reactions) && !length(requested_reactions)) {
      stop(
        paste(
          "The selected genes do not resolve to reactions that were core",
          "targets in the previous LP analysis."
        ),
        call. = FALSE
      )
    }
    gene_source <- lapply(gene_reactions, function(reaction) {
      sort(unique(gpr$gene[
        gpr$reaction_id == reaction & gpr$gene %in% requested_genes
      ]))
    })
    names(gene_source) <- gene_reactions
  }

  reactions <- union(requested_reactions, gene_reactions)
  source <- vapply(reactions, function(reaction) {
    by_id <- reaction %in% requested_reactions
    by_gene <- reaction %in% gene_reactions
    if (by_id && by_gene) {
      "previous_core_reaction_id+gene"
    } else if (by_id) {
      "previous_core_reaction_id"
    } else if (identical(gene_match, "complete_gpr")) {
      "previous_core_gene_complete_gpr"
    } else {
      "previous_core_gene_any_direct"
    }
  }, character(1))
  mapped_genes <- vapply(reactions, function(reaction) {
    genes <- gene_source[[reaction]]
    if (is.null(genes) || !length(genes)) {
      NA_character_
    } else {
      paste(genes, collapse = ";")
    }
  }, character(1))

  data.frame(
    sample_id = "global",
    module_id = "GLOBAL_UNION",
    gene = mapped_genes,
    reaction_id = reactions,
    is_core = TRUE,
    selection_source = source,
    stringsAsFactors = FALSE
  )
}

.rc_build_target_union_definition <- function(
    gem, global_core_reactions, global_reaction_membership,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L) {
  gene_match <- match.arg(gene_match)
  expansion_mode <- match.arg(expansion_mode)
  required <- c("reaction_id")
  if (!is.data.frame(global_core_reactions) ||
      !all(required %in% colnames(global_core_reactions))) {
    stop("`global_core_reactions` must contain reaction_id.", call. = FALSE)
  }
  if (!is.data.frame(global_reaction_membership) ||
      !all(required %in% colnames(global_reaction_membership))) {
    stop("`global_reaction_membership` must contain reaction_id.",
         call. = FALSE)
  }

  selected <- .rc_target_union_core_rows(
    gem = gem,
    available_core_reactions = global_core_reactions$reaction_id,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
  expanded <- rc_expand_meta_module_reactions(
    gem = gem,
    core_reactions = selected,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )
  targets <- expanded$reaction_membership
  previous_union_ids <- .rc_target_union_normalize_ids(
    global_reaction_membership$reaction_id
  )
  missing_from_union <- setdiff(
    unique(as.character(targets$reaction_id)),
    previous_union_ids
  )
  if (length(missing_from_union)) {
    stop(
      paste(
        "Annotation expansion produced reactions absent from the previously",
        "constructed global union GEM:"
      ),
      paste(utils::head(missing_from_union, 10L), collapse = ", "),
      ". Rebuild the original union with matching expansion settings.",
      call. = FALSE
    )
  }

  union_match <- match(
    as.character(targets$reaction_id),
    as.character(global_reaction_membership$reaction_id)
  )
  targets$selected_core_anchor <-
    as.character(targets$reaction_id) %in% selected$reaction_id
  targets$score_target <- TRUE
  targets$target_role <- ifelse(
    targets$selected_core_anchor,
    "selected_previous_core",
    as.character(targets$inclusion_stage)
  )
  targets$previous_union_is_core <- if (
    "is_core" %in% colnames(global_reaction_membership)
  ) {
    global_reaction_membership$is_core[union_match] %in% TRUE
  } else {
    as.character(targets$reaction_id) %in%
      as.character(global_core_reactions$reaction_id)
  }
  targets$previous_union_inclusion_stage <- if (
    "inclusion_stage" %in% colnames(global_reaction_membership)
  ) {
    as.character(global_reaction_membership$inclusion_stage[union_match])
  } else {
    NA_character_
  }
  rownames(selected) <- NULL
  rownames(targets) <- NULL

  summary <- expanded$summary
  summary$n_selected_previous_core <- nrow(selected)
  summary$n_expanded_score_targets <- length(unique(targets$reaction_id))
  summary$n_previous_union_reactions <- length(previous_union_ids)
  summary$gene_match <- gene_match
  summary$expansion_mode <- expansion_mode
  summary$scoring_policy <-
    "selected_core_plus_annotation_related_reactions_all_scored"
  summary$model_policy <- "reuse_previous_global_union_gem_without_rebuilding"

  list(
    selected_core_reactions = selected,
    expanded_scoring_targets = targets,
    previous_union_membership = global_reaction_membership,
    summary = summary,
    crossref_maps = expanded$crossref_maps,
    params = list(
      gene_match = gene_match,
      expansion_mode = expansion_mode,
      max_iterations = as.integer(max_iterations),
      selected_core_reactions = unique(as.character(selected$reaction_id)),
      score_targets = unique(as.character(targets$reaction_id))
    )
  )
}

.rc_build_target_union_model_cache <- function(
    layer2, target_reactions,
    target_direction = c("both", "forward", "reverse")) {
  target_direction <- match.arg(target_direction)
  if (!is.list(layer2) ||
      !identical(as.character(layer2$model_mode), "meta_module_gem")) {
    stop(
      paste(
        "`layer2` must be the completed original core LP result produced with",
        "`model_mode = \"meta_module_gem\"`."
      ),
      call. = FALSE
    )
  }
  summary <- layer2$model_cache_summary
  required <- c("medium_scenario", "file")
  if (!is.data.frame(summary) || !all(required %in% colnames(summary)) ||
      !nrow(summary)) {
    stop(
      "`layer2$model_cache_summary` does not identify reusable union GEM files.",
      call. = FALSE
    )
  }
  summary$medium_scenario <- as.character(summary$medium_scenario)
  summary$file <- as.character(summary$file)
  summary <- unique(summary[
    !is.na(summary$medium_scenario) & nzchar(summary$medium_scenario) &
      !is.na(summary$file) & nzchar(summary$file),
    , drop = FALSE
  ])
  if (!nrow(summary)) {
    stop("No reusable union GEM files remain after cache validation.",
         call. = FALSE)
  }
  scenario_files <- split(summary$file, summary$medium_scenario)
  ambiguous <- names(scenario_files)[vapply(
    scenario_files, function(x) length(unique(x)) != 1L, logical(1)
  )]
  if (length(ambiguous)) {
    stop(
      "Each medium scenario must resolve to one previous union GEM file: ",
      paste(ambiguous, collapse = ", "),
      call. = FALSE
    )
  }
  summary <- summary[
    !duplicated(summary$medium_scenario),
    , drop = FALSE
  ]
  missing_files <- summary$file[!file.exists(summary$file)]
  if (length(missing_files)) {
    stop(
      "Previous union GEM cache files are unavailable: ",
      paste(utils::head(missing_files, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  cache <- list()
  diagnostics <- list()
  for (i in seq_len(nrow(summary))) {
    scenario <- summary$medium_scenario[[i]]
    file <- summary$file[[i]]
    model <- readRDS(file)
    validated <- rc_validate_gem(model)
    missing_targets <- setdiff(target_reactions, validated$reactions)
    if (length(missing_targets)) {
      stop(
        "Scoring targets are absent from the previous union GEM for `",
        scenario, "`: ",
        paste(utils::head(missing_targets, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    directions <- rc_prepare_directional_targets(
      model,
      target_reactions = target_reactions,
      target_direction = target_direction
    )
    directions$medium_scenario <- scenario
    diagnostics[[scenario]] <- directions
    allowed <- directions[
      directions$target_direction %in% c("forward", "reverse"),
      , drop = FALSE
    ]
    for (j in seq_len(nrow(allowed))) {
      reaction <- as.character(allowed$reaction_id[[j]])
      direction <- as.character(allowed$target_direction[[j]])
      key <- paste0(
        "reaction=", utils::URLencode(reaction, reserved = TRUE),
        "::direction=", direction,
        "::medium=", utils::URLencode(scenario, reserved = TRUE)
      )
      cache[[key]] <- list(
        sample_id = "global",
        module_id = "GLOBAL_UNION",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        build_strategy = "reuse_previous_global_union_gem"
      )
    }
  }
  if (!length(cache)) {
    stop("No selected or annotation-related reaction direction can be scored.",
         call. = FALSE)
  }
  summary$build_strategy <- "reuse_previous_global_union_gem"
  summary$reused_without_rebuilding <- TRUE
  attr(cache, "summary") <- summary
  attr(cache, "direction_diagnostics") <- .rc_bind_frames_fill(diagnostics)
  cache
}

.rc_target_union_no_constraint_medium <- function(scenario) {
  data.frame(
    medium_scenario_id = as.character(scenario),
    exchange_reaction_id = NA_character_,
    lb = NA_real_,
    ub = NA_real_,
    available = FALSE,
    .no_constraints = TRUE,
    stringsAsFactors = FALSE
  )
}

.rc_bind_target_union_results <- function(
    results, model_cache, original_medium_scenarios,
    omega, target_direction, solver, time_limit, flux_threshold) {
  if (!length(results)) {
    stop("No second-pass union-GEM results were produced.", call. = FALSE)
  }
  bind_matrix <- function(field) {
    values <- lapply(results, `[[`, field)
    out <- do.call(rbind, values)
    if (is.null(rownames(out)) || anyDuplicated(rownames(out))) {
      stop("Second-pass target row IDs are missing or duplicated.",
           call. = FALSE)
    }
    out
  }
  score <- bind_matrix("score")
  penalty <- bind_matrix("penalty")
  vmax <- bind_matrix("vmax")
  feasible <- bind_matrix("feasible")
  evaluated <- bind_matrix("evaluated")

  target_direction_table <- unique(.rc_bind_frames_fill(Map(
    function(result, scenario) {
      directions <- result$target_direction %||% data.frame()
      if (nrow(directions)) directions$medium_scenario <- scenario
      directions
    },
    results,
    names(results)
  )))
  direction_diagnostics <- .rc_bind_frames_fill(Map(
    function(result, scenario) {
      diagnostics <- result$direction_diagnostics %||% data.frame()
      if (nrow(diagnostics)) diagnostics$medium_scenario <- scenario
      diagnostics
    },
    results,
    names(results)
  ))
  lp_diagnostics <- .rc_bind_frames_fill(lapply(
    results, `[[`, "lp_diagnostics"
  ))

  source_summary <- attr(model_cache, "summary")
  source_files <- stats::setNames(
    as.character(source_summary$file),
    as.character(source_summary$medium_scenario)
  )
  model_cache_summary <- .rc_bind_frames_fill(Map(
    function(result, scenario) {
      summary <- result$model_cache_summary %||% data.frame()
      if (nrow(summary)) {
        summary$source_union_model_file <- unname(source_files[[scenario]])
        summary$reused_previous_union_gem <- TRUE
      }
      summary
    },
    results,
    names(results)
  ))
  original_model_diagnostics <- .rc_bind_frames_fill(lapply(
    seq_len(nrow(source_summary)),
    function(i) {
      file <- as.character(source_summary$file[[i]])
      scenario <- as.character(source_summary$medium_scenario[[i]])
      model <- readRDS(file)
      diagnostics <- model$closure_diagnostics %||% data.frame()
      if (nrow(diagnostics)) {
        diagnostics$medium_scenario <- scenario
        diagnostics$source_union_model_file <- file
        diagnostics$diagnostic_origin <- "original_core_union_build"
      }
      diagnostics
    }
  ))
  components_by_medium <- lapply(
    results,
    `[[`,
    "penalty_components"
  )
  manifests <- .rc_bind_frames_fill(Map(
    function(result, scenario) {
      manifest <- result$model_file_manifest %||% data.frame()
      if (nrow(manifest)) manifest$medium_scenario <- scenario
      manifest
    },
    results,
    names(results)
  ))

  score_semantics <- unique(vapply(
    results,
    function(result) as.character(
      result$score_semantics %||%
        "within_target_relative_penalty_rank_not_probability"
    ),
    character(1)
  ))
  if (length(score_semantics) != 1L) {
    stop("Second-pass score semantics differ across medium scenarios.",
         call. = FALSE)
  }
  noninformative_values <- lapply(results, function(result) {
    value <- result$noninformative_target
    if (is.null(value)) {
      value <- stats::setNames(
        rep(FALSE, nrow(result$score)),
        rownames(result$score)
      )
    }
    value
  })
  noninformative <- unlist(
    unname(noninformative_values),
    use.names = TRUE
  )
  attr(score, "score_semantics") <- score_semantics[[1L]]
  attr(score, "noninformative_target") <- noninformative

  answer <- list(
    score = score,
    penalty = penalty,
    vmax = vmax,
    feasible = feasible,
    evaluated = evaluated,
    target_direction = target_direction_table,
    direction_diagnostics = direction_diagnostics,
    medium_scenarios = original_medium_scenarios,
    model_mode = "reused_global_union_gem",
    model_cache_summary = model_cache_summary,
    model_diagnostics = original_model_diagnostics,
    lp_diagnostics = lp_diagnostics,
    penalty_components = if (length(components_by_medium) == 1L) {
      components_by_medium[[1L]]
    } else {
      components_by_medium
    },
    penalty_components_by_medium = components_by_medium,
    evidence_policy = results[[1L]]$evidence_policy,
    evidence_policy_detail = results[[1L]]$evidence_policy_detail,
    unit_meta = results[[1L]]$unit_meta,
    params = list(
      unit = "metacell",
      omega = omega,
      target_direction = target_direction,
      shared_gem = TRUE,
      shared_gem_scope = "previous_global_union_by_medium",
      reused_without_structural_reconstruction = TRUE,
      second_pass_engine = "canonical_full_gem_microcompass",
      parallel_task = "reused_union_model_by_metacell",
      flux_threshold = flux_threshold,
      solver = solver,
      time_limit = time_limit
    ),
    method = paste(
      "canonical microCOMPASS directional LP on previously constructed",
      "global union GEMs"
    )
  )
  if (nrow(manifests)) answer$model_file_manifest <- manifests
  answer$relative_penalty_rank <- answer$score
  answer$score_semantics <- score_semantics[[1L]]
  answer$noninformative_target <- noninformative
  answer$primary_output <- "penalty"
  answer$primary_output_semantics <-
    "minimum evidence-discordance penalty; lower means stronger support"
  answer
}

.rc_score_target_union_cache <- function(
    layer1, gem, model_cache, medium_scenarios,
    condition_col, sample_col, celltype_col,
    omega = 0.95,
    target_direction = c("both", "forward", "reverse"),
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
    cache_dir = tempfile("RegCompassR_target_union_"),
    parallel = TRUE, BPPARAM = NULL) {
  target_direction <- match.arg(target_direction)
  solver <- match.arg(solver)
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must be one finite value in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(time_limit) || length(time_limit) != 1L ||
      !is.finite(time_limit) || time_limit <= 0) {
    stop("`time_limit` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(flux_threshold) || length(flux_threshold) != 1L ||
      !is.finite(flux_threshold) || flux_threshold < 0) {
    stop("`flux_threshold` must be one finite non-negative number.",
         call. = FALSE)
  }
  if (!is.character(cache_dir) || length(cache_dir) != 1L ||
      is.na(cache_dir) || !nzchar(cache_dir)) {
    stop("`cache_dir` must be one non-empty path.", call. = FALSE)
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  source_summary <- attr(model_cache, "summary")
  results <- list()
  safe <- function(value) {
    paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))), collapse = "")
  }
  for (i in seq_len(nrow(source_summary))) {
    scenario <- as.character(source_summary$medium_scenario[[i]])
    source_file <- as.character(source_summary$file[[i]])
    source_model <- readRDS(source_file)
    scenario_entries <- model_cache[vapply(
      model_cache,
      function(entry) identical(
        as.character(entry$medium_scenario),
        scenario
      ),
      logical(1)
    )]
    if (!length(scenario_entries)) next
    target_reactions <- unique(vapply(
      scenario_entries,
      `[[`,
      character(1),
      "reaction_id"
    ))
    scenario_cache_dir <- file.path(
      cache_dir,
      paste0("medium_", safe(scenario))
    )
    results[[scenario]] <- withCallingHandlers(
      rc_run_microcompass(
        layer1 = layer1,
        gem = source_model,
        target_reactions = target_reactions,
        medium_scenarios = .rc_target_union_no_constraint_medium(scenario),
        mode = "full_gem",
        unit = "metacell",
        condition_col = condition_col,
        sample_col = sample_col,
        celltype_col = celltype_col,
        model_params = list(cache_dir = scenario_cache_dir),
        omega = omega,
        target_direction = target_direction,
        parallel = parallel,
        solver = solver,
        time_limit = time_limit,
        flux_threshold = flux_threshold,
        BPPARAM = BPPARAM
      ),
      warning = function(w) {
        if (grepl(
          "Metacell-level scores are descriptive pseudo-observations",
          conditionMessage(w),
          fixed = TRUE
        )) invokeRestart("muffleWarning")
      }
    )
  }
  .rc_bind_target_union_results(
    results = results,
    model_cache = model_cache,
    original_medium_scenarios = medium_scenarios,
    omega = omega,
    target_direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold
  )
}

#' Re-score annotation-related reactions in a previous global union GEM
#'
#' This optional post-Layer-2 step starts from reactions that were already core
#' targets in the original LP analysis. It resolves those cores from reaction IDs
#' and/or GPR genes, expands them through same-subsystem, shared KEGG/Reactome
#' reaction, and shared master-Rhea mappings, and scores every expanded reaction
#' as an LP target. The structural model is loaded from the previous Layer-2
#' global union-GEM cache and is not reconstructed.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer2 Completed original core LP output from
#'   [rc_regcompass_step_layer2()] using `model_mode = "meta_module_gem"`.
#' @param gem The same validated GEM used by the original analysis.
#' @param outdir Persistent output directory.
#' @param core_reaction_ids One or more reactions that were core targets in the
#'   previous LP analysis.
#' @param core_genes Genes used to resolve previous core reactions through GPR
#'   rules.
#' @param gene_match `"complete_gpr"` requires one complete GPR AND group;
#'   `"any_direct"` permits an intentional partial-complex match.
#' @param subsystem_table Optional external reaction-to-subsystem table.
#' @param expansion_mode Annotation expansion mode passed to
#'   [rc_expand_meta_module_reactions()].
#' @param max_iterations Maximum fixed-point annotation-expansion iterations.
#' @param layer2_args Optional `omega`, `target_direction`, `solver`,
#'   `time_limit`, and `flux_threshold` overrides for the second scoring pass.
#' @param parallel Whether LP tasks may run in parallel.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @return A `regcompass_target_union_step` containing selected previous cores,
#'   expanded LP targets, the reused union membership, and the second-pass
#'   microCOMPASS result.
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L,
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {
  gene_match <- match.arg(gene_match)
  expansion_mode <- match.arg(expansion_mode)
  if (!is.list(layer1) || is.null(layer1$reaction_expression) ||
      is.null(layer1$unit_meta)) {
    stop("`layer1` must contain reaction_expression and unit_meta.",
         call. = FALSE)
  }
  if (!inherits(meta_modules, "regcompass_meta_module_step")) {
    stop(
      "`meta_modules` must be the output of `rc_regcompass_step_meta_modules()`.",
      call. = FALSE
    )
  }
  if (!is.list(layer2) || is.null(layer2$penalty)) {
    stop(
      "`layer2` must be the completed original core LP result.",
      call. = FALSE
    )
  }
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c(
    "omega", "target_direction", "solver",
    "time_limit", "flux_threshold"
  )
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop(
      "Unsupported `layer2_args` for union-GEM re-scoring: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  global <- meta_modules$global_modules
  if (!is.list(global) ||
      !is.data.frame(global$global_core_reactions) ||
      !is.data.frame(global$global_reaction_membership)) {
    stop("The previous global meta-module union is unavailable.",
         call. = FALSE)
  }
  definition <- .rc_build_target_union_definition(
    gem = gem,
    global_core_reactions = global$global_core_reactions,
    global_reaction_membership = global$global_reaction_membership,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )

  target_direction <- match.arg(
    as.character(
      layer2_args$target_direction %||%
        layer2$params$target_direction %||%
        "both"
    ),
    c("both", "forward", "reverse")
  )
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  omega <- layer2_args$omega %||% layer2$params$omega %||% 0.95
  time_limit <- layer2_args$time_limit %||% 60
  flux_threshold <- layer2_args$flux_threshold %||% 1e-8

  model_cache <- .rc_build_target_union_model_cache(
    layer2 = layer2,
    target_reactions = definition$params$score_targets,
    target_direction = target_direction
  )
  workflow <- meta_modules$workflow_params
  scored <- .rc_score_target_union_cache(
    layer1 = layer1,
    gem = gem,
    model_cache = model_cache,
    medium_scenarios = layer2$medium_scenarios,
    condition_col = workflow$condition_col,
    sample_col = workflow$sample_col,
    celltype_col = workflow$celltype_col,
    omega = omega,
    target_direction = target_direction,
    solver = solver,
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    cache_dir = file.path(outdir, "model_cache"),
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  scored$params$target_scope <-
    "selected_previous_core_plus_annotation_related_reactions"
  scored$params$n_selected_previous_core <-
    nrow(definition$selected_core_reactions)
  scored$params$n_expanded_score_targets <-
    length(definition$params$score_targets)
  scored$params$annotation_expansion <- c(
    "same_subsystem",
    "shared_kegg_or_reactome_reaction",
    "shared_master_rhea_reaction"
  )

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_mm_write_tsv_gz(
    definition$selected_core_reactions,
    file.path(outdir, "selected_previous_core_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$expanded_scoring_targets,
    file.path(outdir, "expanded_scoring_targets.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$previous_union_membership,
    file.path(outdir, "reused_global_union_membership.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$summary,
    file.path(outdir, "target_union_summary.tsv.gz")
  )
  rc_export_microcompass(scored, file.path(outdir, "scores"))
  answer <- c(definition, list(microcompass = scored))
  class(answer) <- c("regcompass_target_union_step", "list")
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
