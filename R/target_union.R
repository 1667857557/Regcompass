# Score direct database-linked non-core reactions in final union GEMs.

.rc_target_union_normalize_ids <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

.rc_target_union_direct_crossref_relations <- function(
    gem, selected_core_reactions) {
  maps <- rc_reaction_crossref_maps(gem)
  specifications <- list(
    list(
      map = .rc_clean_meta_module_map(maps$kegg, "kegg_id"),
      id_col = "kegg_id",
      expansion_type = "shared_kegg_reaction",
      prefix = "KEGG:"
    ),
    list(
      map = .rc_clean_meta_module_map(maps$reactome, "reactome_id"),
      id_col = "reactome_id",
      expansion_type = "shared_reactome_reaction",
      prefix = "REACTOME:"
    ),
    list(
      map = .rc_clean_meta_module_map(
        maps$rhea_master, "rhea_master_id"
      ),
      id_col = "rhea_master_id",
      expansion_type = "shared_master_rhea_reaction",
      prefix = "RHEA_MASTER:"
    )
  )
  anchors <- .rc_target_union_normalize_ids(
    selected_core_reactions$reaction_id
  )
  output <- list()
  output_index <- 0L
  for (anchor in anchors) {
    for (specification in specifications) {
      map <- specification$map
      id_col <- specification$id_col
      if (!is.data.frame(map) || !nrow(map)) next
      anchor_ids <- unique(as.character(map[[id_col]][
        map$reaction_id == anchor
      ]))
      anchor_ids <- anchor_ids[
        !is.na(anchor_ids) & nzchar(trimws(anchor_ids))
      ]
      if (!length(anchor_ids)) next
      reactions <- unique(as.character(map$reaction_id[
        map[[id_col]] %in% anchor_ids
      ]))
      reactions <- setdiff(reactions, anchor)
      for (reaction in reactions) {
        shared_ids <- intersect(
          anchor_ids,
          unique(as.character(map[[id_col]][map$reaction_id == reaction]))
        )
        shared_ids <- sort(shared_ids[
          !is.na(shared_ids) & nzchar(trimws(shared_ids))
        ])
        if (!length(shared_ids)) next
        output_index <- output_index + 1L
        output[[output_index]] <- data.frame(
          anchor_core_reaction_id = anchor,
          reaction_id = reaction,
          expansion_type = specification$expansion_type,
          source_annotation = paste0(
            specification$prefix, paste(shared_ids, collapse = ";")
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(output)) {
    return(data.frame(
      anchor_core_reaction_id = character(),
      reaction_id = character(),
      expansion_type = character(),
      source_annotation = character(),
      stringsAsFactors = FALSE
    ))
  }
  answer <- unique(do.call(rbind, output))
  answer <- answer[order(
    answer$anchor_core_reaction_id,
    answer$reaction_id,
    answer$expansion_type,
    answer$source_annotation
  ), , drop = FALSE]
  rownames(answer) <- NULL
  answer
}

.rc_target_union_aggregate_targets <- function(catalogue) {
  rows <- split(seq_len(nrow(catalogue)), catalogue$reaction_id)
  answer <- do.call(rbind, lapply(rows, function(index) {
    one <- catalogue[index, , drop = FALSE]
    stages <- sort(unique(as.character(
      one$merged_catalogue_inclusion_stage[
        !is.na(one$merged_catalogue_inclusion_stage) &
          nzchar(one$merged_catalogue_inclusion_stage)
      ]
    )))
    data.frame(
      catalogue_id = "MERGED_META_MODULES",
      reaction_id = as.character(one$reaction_id[[1L]]),
      anchor_core_reaction_ids = paste(
        sort(unique(as.character(one$anchor_core_reaction_id))),
        collapse = ";"
      ),
      expansion_types = paste(
        sort(unique(as.character(one$expansion_type))),
        collapse = ";"
      ),
      source_annotations = paste(
        sort(unique(as.character(one$source_annotation))),
        collapse = ";"
      ),
      merged_catalogue_is_core = FALSE,
      merged_catalogue_inclusion_stage = if (length(stages)) {
        paste(stages, collapse = ";")
      } else {
        NA_character_
      },
      score_target = TRUE,
      target_role = "direct_database_crossref_noncore",
      lp_exclusion_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(answer) <- NULL
  answer
}

.rc_target_union_model_summary <- function(layer2) {
  if (!inherits(layer2, "regcompass_layer2_step") ||
      !identical(as.character(layer2$model_mode), "meta_module_gem")) {
    stop(
      paste(
        "`layer2` must be the completed core LP stage with",
        "`model_mode = \"meta_module_gem\"`."
      ),
      call. = FALSE
    )
  }
  summary <- layer2$model_cache_summary
  required <- c(
    "medium_scenario", "file", "file_checksum",
    "build_strategy", "completion_stage"
  )
  if (!is.data.frame(summary) || !all(required %in% colnames(summary)) ||
      !nrow(summary)) {
    stop(
      "`layer2$model_cache_summary` does not identify final union GEM files.",
      call. = FALSE
    )
  }
  summary$medium_scenario <- trimws(as.character(summary$medium_scenario))
  summary$file <- as.character(summary$file)
  summary$file_checksum <- as.character(summary$file_checksum)
  summary <- unique(summary[
    !is.na(summary$medium_scenario) & nzchar(summary$medium_scenario) &
      !is.na(summary$file) & nzchar(summary$file),
    , drop = FALSE
  ])
  scenario_files <- split(summary$file, summary$medium_scenario)
  ambiguous <- names(scenario_files)[vapply(
    scenario_files, function(x) length(unique(x)) != 1L, logical(1)
  )]
  if (length(ambiguous)) {
    stop(
      "Each medium scenario must resolve to one final union GEM file: ",
      paste(ambiguous, collapse = ", "),
      call. = FALSE
    )
  }
  summary <- summary[!duplicated(summary$medium_scenario), , drop = FALSE]
  for (i in seq_len(nrow(summary))) {
    if (!identical(
      as.character(summary$build_strategy[[i]]),
      "medium_specific_union_gem"
    ) || !identical(
      as.character(summary$completion_stage[[i]]),
      "single_global_fastcore_after_meta_module_merge"
    )) {
      stop(
        "Second-pass scoring requires Stage 5 final medium-specific union GEMs.",
        call. = FALSE
      )
    }
    .rc_read_cached_union_gem(
      file = summary$file[[i]],
      medium_scenario = summary$medium_scenario[[i]],
      expected_checksum = summary$file_checksum[[i]]
    )
  }
  summary
}

.rc_target_union_cached_reaction_ids <- function(layer2) {
  summary <- .rc_target_union_model_summary(layer2)
  reaction_sets <- lapply(seq_len(nrow(summary)), function(i) {
    model <- .rc_read_cached_union_gem(
      file = summary$file[[i]],
      medium_scenario = summary$medium_scenario[[i]],
      expected_checksum = summary$file_checksum[[i]]
    )
    rc_validate_gem(model)$reactions
  })
  available <- Reduce(intersect, reaction_sets)
  available <- .rc_target_union_normalize_ids(available)
  if (!length(available)) {
    stop("The final medium-specific union GEMs share no reactions.",
         call. = FALSE)
  }
  attr(available, "model_summary") <- summary
  attr(available, "reaction_counts") <- vapply(
    reaction_sets, length, integer(1)
  )
  available
}

.rc_build_target_union_definition <- function(
    gem, merged_core_reactions, merged_reaction_membership,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    cached_reaction_ids = NULL) {
  gene_match <- match.arg(gene_match)
  if (!is.data.frame(merged_core_reactions) ||
      !"reaction_id" %in% colnames(merged_core_reactions)) {
    stop("`merged_core_reactions` must contain reaction_id.",
         call. = FALSE)
  }
  if (!is.data.frame(merged_reaction_membership) ||
      !"reaction_id" %in% colnames(merged_reaction_membership)) {
    stop("`merged_reaction_membership` must contain reaction_id.",
         call. = FALSE)
  }
  selected <- .rc_target_union_core_rows(
    gem = gem,
    available_core_reactions = merged_core_reactions$reaction_id,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
  catalogue <- .rc_target_union_direct_crossref_relations(gem, selected)
  if (!nrow(catalogue)) {
    stop(
      paste(
        "The selected core reactions have no directly linked KEGG, Reactome,",
        "or master-Rhea reactions."
      ),
      call. = FALSE
    )
  }
  membership_ids <- .rc_target_union_normalize_ids(
    merged_reaction_membership$reaction_id
  )
  available_ids <- .rc_target_union_normalize_ids(cached_reaction_ids)
  if (!length(available_ids)) {
    stop("No reusable reactions were found in the final union GEM cache.",
         call. = FALSE)
  }
  catalogue_match <- match(
    as.character(catalogue$reaction_id),
    as.character(merged_reaction_membership$reaction_id)
  )
  catalogue$available_in_all_cached_union_gems <-
    catalogue$reaction_id %in% available_ids
  catalogue$present_in_merged_catalogue <- !is.na(catalogue_match)
  catalogue$merged_catalogue_is_core <- catalogue$reaction_id %in%
    as.character(merged_core_reactions$reaction_id)
  catalogue$merged_catalogue_inclusion_stage <- if (
    "inclusion_stage" %in% colnames(merged_reaction_membership)
  ) {
    as.character(merged_reaction_membership$inclusion_stage[catalogue_match])
  } else {
    rep(NA_character_, nrow(catalogue))
  }
  support_only <- catalogue$available_in_all_cached_union_gems &
    !catalogue$present_in_merged_catalogue
  catalogue$merged_catalogue_inclusion_stage[support_only] <-
    "global_fastcore_support_not_in_merged_catalogue"
  catalogue$score_target <- !catalogue$merged_catalogue_is_core &
    catalogue$available_in_all_cached_union_gems
  catalogue$target_role <- ifelse(
    catalogue$merged_catalogue_is_core,
    "merged_core_not_rescored",
    ifelse(
      catalogue$available_in_all_cached_union_gems,
      "direct_database_crossref_noncore",
      "direct_database_crossref_absent_from_cached_union_gem"
    )
  )
  catalogue$lp_exclusion_reason <- ifelse(
    catalogue$merged_catalogue_is_core,
    "already_scored_in_original_layer2",
    ifelse(
      catalogue$available_in_all_cached_union_gems,
      NA_character_,
      "absent_from_one_or_more_cached_union_gems"
    )
  )
  target_relations <- catalogue[catalogue$score_target, , drop = FALSE]
  if (!nrow(target_relations)) {
    stop(
      paste(
        "No directly linked non-core reaction is present in every final",
        "medium-specific union GEM. Original core reactions are not recomputed."
      ),
      call. = FALSE
    )
  }
  targets <- .rc_target_union_aggregate_targets(target_relations)
  rownames(selected) <- NULL
  rownames(catalogue) <- NULL
  summary <- data.frame(
    n_selected_core = nrow(selected),
    n_direct_crossref_relations = nrow(catalogue),
    n_direct_crossref_reactions = length(unique(catalogue$reaction_id)),
    n_merged_core_reactions_not_rescored = length(unique(
      catalogue$reaction_id[catalogue$merged_catalogue_is_core]
    )),
    n_cached_union_unavailable_reactions = length(unique(
      catalogue$reaction_id[!catalogue$available_in_all_cached_union_gems]
    )),
    n_expanded_score_targets = nrow(targets),
    n_merged_catalogue_reactions = length(membership_ids),
    n_reactions_shared_by_cached_union_gems = length(available_ids),
    gene_match = gene_match,
    expansion_policy =
      "direct_from_selected_core_via_kegg_reactome_master_rhea_only",
    scoring_policy =
      "direct_noncore_reactions_present_in_all_final_union_gems",
    model_policy = "reuse_exact_final_medium_specific_union_gems",
    stringsAsFactors = FALSE
  )
  list(
    selected_core_reactions = selected,
    expanded_reaction_catalog = catalogue,
    expanded_scoring_targets = targets,
    merged_catalogue_membership = merged_reaction_membership,
    summary = summary,
    params = list(
      gene_match = gene_match,
      selected_core_reactions = unique(as.character(selected$reaction_id)),
      merged_core_reactions_not_rescored = unique(as.character(
        catalogue$reaction_id[catalogue$merged_catalogue_is_core]
      )),
      cached_union_unavailable_reactions = unique(as.character(
        catalogue$reaction_id[!catalogue$available_in_all_cached_union_gems]
      )),
      score_targets = unique(as.character(targets$reaction_id)),
      expansion_policy =
        "direct_from_selected_core_via_kegg_reactome_master_rhea_only"
    )
  )
}

.rc_build_target_union_model_cache <- function(
    layer2, target_reactions,
    target_direction = c("both", "forward", "reverse")) {
  target_direction <- match.arg(target_direction)
  summary <- .rc_target_union_model_summary(layer2)
  cache <- list()
  diagnostics <- list()
  fingerprints <- character(nrow(summary))
  for (i in seq_len(nrow(summary))) {
    scenario <- summary$medium_scenario[[i]]
    file <- summary$file[[i]]
    checksum <- summary$file_checksum[[i]]
    model <- .rc_read_cached_union_gem(file, scenario, checksum)
    validated <- rc_validate_gem(model)
    fingerprints[[i]] <- .rc_full_gem_cache_fingerprint(model)
    missing_targets <- setdiff(target_reactions, validated$reactions)
    if (length(missing_targets)) {
      stop(
        "Scoring targets are absent from the final union GEM for `",
        scenario, "`: ",
        paste(utils::head(missing_targets, 10L), collapse = ", "),
        call. = FALSE
      )
    }
    directions <- rc_prepare_directional_targets(
      model, target_reactions, target_direction
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
        module_id = "MEDIUM_UNION_GEM",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        file_checksum = checksum,
        build_strategy = "reuse_exact_final_medium_specific_union_gem"
      )
    }
  }
  if (!length(cache)) {
    stop("No direct database-linked non-core reaction direction can be scored.",
         call. = FALSE)
  }
  summary$source_model_fingerprint <- fingerprints
  summary$reused_without_rebuilding <- TRUE
  summary$second_pass_build_strategy <-
    "reuse_exact_final_medium_specific_union_gem"
  attr(cache, "summary") <- summary
  attr(cache, "direction_diagnostics") <- .rc_bind_frames_fill(diagnostics)
  cache
}

.rc_score_existing_union_cache <- function(
    layer1, gem, model_cache,
    condition_col, sample_col, celltype_col,
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    flux_threshold = 1e-8,
    parallel = TRUE, BPPARAM = NULL) {
  solver <- match.arg(solver)
  .rc_require_lp_solver(solver)
  if (!is.numeric(omega) || length(omega) != 1L ||
      !is.finite(omega) || omega <= 0 || omega > 1) {
    stop("`omega` must be one finite value in (0, 1].", call. = FALSE)
  }
  matrices <- rc_layer2_unit_matrices(
    layer1, "metacell", sample_col, celltype_col, condition_col
  )
  row_ids <- names(model_cache)
  units <- colnames(matrices$reaction_expression)
  model_files <- vapply(model_cache, `[[`, character(1), "file")
  unique_files <- unique(model_files)
  representative <- vapply(unique_files, function(file) {
    row_ids[match(file, model_files)]
  }, character(1))
  all_reactions <- unique(unlist(lapply(representative, function(row_id) {
    entry <- model_cache[[row_id]]
    model <- .rc_read_cached_union_gem(
      entry$file, entry$medium_scenario, entry$file_checksum
    )
    colnames(model$S)
  }), use.names = FALSE))
  gem <- rc_annotate_reaction_roles(gem)
  penalties <- rc_compute_multiome_penalty(
    rc_align_reaction_expression(
      matrices$reaction_expression, all_reactions, NA_real_
    ),
    reaction_roles = gem$reaction_roles
  )
  penalty <- vmax <- matrix(
    NA_real_, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )
  feasible <- evaluated <- matrix(
    FALSE, length(row_ids), length(units),
    dimnames = list(row_ids, units)
  )
  tasks <- expand.grid(
    file = unique_files, unit_id = units, stringsAsFactors = FALSE
  )
  run_one <- function(task) {
    file <- as.character(task$file)
    unit_id <- as.character(task$unit_id)
    selected <- row_ids[model_files == file]
    first <- model_cache[[selected[[1L]]]]
    model <- .rc_read_cached_union_gem(
      first$file, first$medium_scenario, first$file_checksum
    )
    answers <- lapply(selected, function(row_id) {
      entry <- model_cache[[row_id]]
      fit <- rc_compass_two_step_lp_directional(
        S = model$S,
        lb = model$lb,
        ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = penalties$penalty[colnames(model$S), unit_id],
        target_direction = entry$target_direction,
        omega = omega,
        solver = solver,
        flux_threshold = flux_threshold
      )
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = fit$penalty,
        vmax = fit$vmax,
        feasible = isTRUE(fit$feasible),
        diagnostics = data.frame(
          row_id = row_id,
          unit_id = unit_id,
          sample_id = "global",
          module_id = "MEDIUM_UNION_GEM",
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
          objective_value = fit$penalty,
          vmax = fit$vmax,
          source_union_gem_file = file,
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      results = answers,
      diagnostics = do.call(rbind, lapply(answers, `[[`, "diagnostics"))
    )
  }
  grouped <- rc_parallel_lapply(
    split(tasks, seq_len(nrow(tasks))),
    function(task) run_one(task[1L, , drop = FALSE]),
    BPPARAM = if (isTRUE(parallel)) BPPARAM else FALSE
  )
  results <- unlist(lapply(grouped, `[[`, "results"), recursive = FALSE)
  for (result in results) {
    penalty[result$row_id, result$unit_id] <- result$penalty
    vmax[result$row_id, result$unit_id] <- result$vmax
    feasible[result$row_id, result$unit_id] <- result$feasible
    evaluated[result$row_id, result$unit_id] <- TRUE
  }
  score <- rc_compass_score_from_penalty(penalty, feasible)
  summary <- attr(model_cache, "summary")
  model_diagnostics <- .rc_bind_frames_fill(lapply(
    seq_len(nrow(summary)), function(i) {
      model <- .rc_read_cached_union_gem(
        summary$file[[i]], summary$medium_scenario[[i]],
        summary$file_checksum[[i]]
      )
      out <- model$closure_diagnostics %||% data.frame()
      if (nrow(out)) {
        out$medium_scenario <- summary$medium_scenario[[i]]
        out$source_union_gem_file <- summary$file[[i]]
      }
      out
    }
  ))
  directions <- unique(do.call(rbind, lapply(model_cache, function(entry) {
    data.frame(
      reaction_id = entry$reaction_id,
      target_direction = entry$target_direction,
      medium_scenario = entry$medium_scenario,
      stringsAsFactors = FALSE
    )
  })))
  manifest <- data.frame(
    file = summary$file,
    medium_scenario = summary$medium_scenario,
    size_bytes = as.numeric(file.info(summary$file)$size),
    checksum = unname(tools::md5sum(summary$file)),
    stringsAsFactors = FALSE
  )
  answer <- list(
    score = score,
    penalty = penalty,
    vmax = vmax,
    feasible = feasible,
    evaluated = evaluated,
    target_direction = directions,
    direction_diagnostics = attr(model_cache, "direction_diagnostics"),
    model_mode = "reused_medium_specific_union_gem",
    model_cache_summary = summary,
    model_diagnostics = model_diagnostics,
    model_file_manifest = manifest,
    lp_diagnostics = .rc_bind_frames_fill(lapply(grouped, `[[`, "diagnostics")),
    penalty_components = penalties$components,
    evidence_policy = penalties$evidence_policy,
    evidence_policy_detail = penalties$evidence_policy_detail,
    unit_meta = matrices$unit_meta,
    params = list(
      unit = "metacell",
      omega = omega,
      shared_gem = TRUE,
      shared_gem_scope = "cached_medium_specific_union_gem_by_medium",
      structural_model_reused_exactly = TRUE,
      fastcore_rerun = FALSE,
      model_rebuild = FALSE,
      parallel_task = "reused_union_gem_by_metacell",
      flux_threshold = flux_threshold,
      solver = solver,
      scoring_time_limit = "none"
    ),
    method = paste(
      "microCOMPASS directional LP for direct",
      "KEGG/Reactome/master-Rhea-linked non-core reactions on exact cached",
      "final medium-specific union GEMs"
    )
  )
  answer$relative_penalty_rank <- answer$score
  answer$score_semantics <- attr(answer$score, "score_semantics") %||%
    "within_target_relative_penalty_rank_not_probability"
  answer$noninformative_target <- attr(answer$score, "noninformative_target")
  answer$primary_output <- "penalty"
  answer$primary_output_semantics <-
    "minimum evidence-discordance penalty; lower means stronger support"
  class(answer) <- c("regcompass_expanded_layer2_result", "list")
  answer
}

#' Score directly database-linked non-core reactions in final union GEMs
#'
#' Selected original core reactions are mapping anchors only. The function
#' scores directly KEGG-, Reactome-, or master-Rhea-linked non-core reactions
#' by reusing the exact final medium-specific union GEM files created by Stage
#' 5. It does not rebuild a GEM and does not rerun FASTCORE.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer2 Output from [rc_regcompass_step_layer2()] with
#'   `model_mode = "meta_module_gem"`.
#' @param gem The same GEM used for the original run.
#' @param outdir Output directory.
#' @param core_reaction_ids GEM reaction IDs used as direct mapping anchors.
#'   The historical argument name is retained for compatibility; an ID may be
#'   an original Layer 2 core or any other valid reaction in `gem`.
#' @param core_genes Genes used to resolve original core anchors through GPRs.
#'   Gene selection intentionally remains restricted to original cores.
#' @param gene_match Require a complete GPR group or allow any direct gene match.
#' @param layer2_args Optional `omega`, `target_direction`, `solver`, and
#'   `flux_threshold` overrides. Scoring LPs have no time-limit control.
#' @param parallel Whether to parallelize model-by-metacell tasks.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @param progress Whether to display stage progress.
#' @return A `regcompass_target_union_step` containing selected anchors, direct
#'   database-linked non-core targets, final union-GEM provenance, and LP scores.
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("target_union", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  gene_match <- match.arg(gene_match)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  workflow <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = workflow, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = workflow, gem = gem,
    required_mode = "meta_module_gem", argument = "layer2"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  allowed <- c("omega", "target_direction", "solver", "flux_threshold")
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop(
      "Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
      ". Scoring `time_limit` has been removed; only ",
      "`layer2_args$model_params$completion_time_limit` in the original ",
      "union-GEM construction stage is supported.",
      call. = FALSE
    )
  }
  catalogue <- meta_modules$merged_modules
  if (!is.list(catalogue) ||
      !is.data.frame(catalogue$merged_core_reactions) ||
      !is.data.frame(catalogue$merged_reaction_membership)) {
    stop("The merged biological meta-module catalogue is unavailable.",
         call. = FALSE)
  }
  if (!setequal(
    as.character(layer2$source_core_reactions$reaction_id),
    as.character(catalogue$merged_core_reactions$reaction_id)
  )) {
    stop("Layer 2 was not generated from the supplied merged core reactions.",
         call. = FALSE)
  }
  cached_reaction_ids <- .rc_target_union_cached_reaction_ids(layer2)
  definition <- .rc_build_target_union_definition(
    gem = gem,
    merged_core_reactions = catalogue$merged_core_reactions,
    merged_reaction_membership = catalogue$merged_reaction_membership,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match,
    cached_reaction_ids = cached_reaction_ids
  )
  target_direction <- match.arg(
    as.character(layer2_args$target_direction %||%
                   layer2$params$target_direction %||% "both"),
    c("both", "forward", "reverse")
  )
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  omega <- layer2_args$omega %||% layer2$params$omega %||% 0.95
  flux_threshold <- layer2_args$flux_threshold %||% 1e-8
  model_cache <- .rc_build_target_union_model_cache(
    layer2 = layer2,
    target_reactions = definition$params$score_targets,
    target_direction = target_direction
  )
  scored <- .rc_score_existing_union_cache(
    layer1 = layer1,
    gem = gem,
    model_cache = model_cache,
    condition_col = workflow$condition_col,
    sample_col = workflow$sample_col,
    celltype_col = workflow$celltype_col,
    omega = omega,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  scored$workflow_params <- workflow
  scored$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  scored$params$target_direction <- target_direction
  scored$params$target_scope <-
    "direct_kegg_reactome_master_rhea_noncore_only"
  scored$params$n_selected_core <- nrow(definition$selected_core_reactions)
  scored$params$n_merged_core_reactions_not_rescored <-
    length(definition$params$merged_core_reactions_not_rescored)
  scored$params$n_cached_union_unavailable_reactions <-
    length(definition$params$cached_union_unavailable_reactions)
  scored$params$n_expanded_score_targets <-
    length(definition$params$score_targets)

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_mm_write_tsv_gz(
    definition$selected_core_reactions,
    file.path(outdir, "selected_core_reactions.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$expanded_reaction_catalog,
    file.path(outdir, "expanded_reaction_catalog.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$expanded_scoring_targets,
    file.path(outdir, "expanded_scoring_targets.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$merged_catalogue_membership,
    file.path(outdir, "merged_meta_module_catalogue_membership.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    definition$summary,
    file.path(outdir, "target_union_summary.tsv.gz")
  )
  rc_export_microcompass(scored, file.path(outdir, "scores"))
  answer <- c(definition, list(microcompass = scored))
  answer$workflow_params <- workflow
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(answer) <- c("regcompass_target_union_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
