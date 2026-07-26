# Active second-pass target-union contract. It reuses exact Stage 5 union GEMs
# and performs metacell-only scoring without sample placeholders.

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
        module_id = "MEDIUM_UNION_GEM",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
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
    layer1, gem, model_cache, condition_col, celltype_col,
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
    layer1,
    celltype_col = celltype_col,
    condition_col = condition_col
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
          module_id = "MEDIUM_UNION_GEM",
          reaction_id = entry$reaction_id,
          target_direction = entry$target_direction,
          medium_scenario = entry$medium_scenario,
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
      aggregation = "none",
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

#' Score direct database-linked reactions without sample aggregation
#'
#' @inheritParams rc_regcompass_step_target_union
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
      ". Scoring `time_limit` and sample aggregation have been removed.",
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
    celltype_col = workflow$celltype_col,
    omega = omega,
    solver = solver,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  scored$workflow_params <- workflow
  scored$workflow_params$sample_col <- NULL
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
  answer$workflow_params$sample_col <- NULL
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  class(answer) <- c("regcompass_target_union_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
