# Score direct database-linked non-core reactions in an existing union GEM.

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
        "Gene-selected cores require a GEM `gpr_table` containing reaction_id, and_group_id and gene.",
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
        "The selected genes do not resolve to core targets in the previous LP analysis.",
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
    if (is.null(genes) || !length(genes)) NA_character_ else
      paste(genes, collapse = ";")
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

.rc_target_union_direct_crossref_relations <- function(gem, selected_core_reactions) {
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
            specification$prefix,
            paste(shared_ids, collapse = ";")
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

.rc_target_union_aggregate_targets <- function(catalog) {
  rows <- split(seq_len(nrow(catalog)), catalog$reaction_id)
  answer <- do.call(rbind, lapply(rows, function(index) {
    one <- catalog[index, , drop = FALSE]
    data.frame(
      sample_id = "global",
      module_id = "GLOBAL_UNION",
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
      previous_union_is_core = FALSE,
      previous_union_inclusion_stage = paste(
        sort(unique(as.character(
          one$previous_union_inclusion_stage[
            !is.na(one$previous_union_inclusion_stage) &
              nzchar(one$previous_union_inclusion_stage)
          ]
        ))),
        collapse = ";"
      ),
      score_target = TRUE,
      target_role = "direct_database_crossref_noncore",
      lp_exclusion_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(answer) <- NULL
  answer
}

.rc_build_target_union_definition <- function(
    gem, global_core_reactions, global_reaction_membership,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct")) {
  gene_match <- match.arg(gene_match)
  if (!is.data.frame(global_core_reactions) ||
      !"reaction_id" %in% colnames(global_core_reactions)) {
    stop("`global_core_reactions` must contain reaction_id.", call. = FALSE)
  }
  if (!is.data.frame(global_reaction_membership) ||
      !"reaction_id" %in% colnames(global_reaction_membership)) {
    stop("`global_reaction_membership` must contain reaction_id.", call. = FALSE)
  }
  selected <- .rc_target_union_core_rows(
    gem = gem,
    available_core_reactions = global_core_reactions$reaction_id,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
  )
  catalog <- .rc_target_union_direct_crossref_relations(gem, selected)
  if (!nrow(catalog)) {
    stop(
      "The selected core reactions have no directly linked KEGG, Reactome, or master-Rhea reactions.",
      call. = FALSE
    )
  }
  previous_union_ids <- .rc_target_union_normalize_ids(
    global_reaction_membership$reaction_id
  )
  missing_from_union <- setdiff(
    unique(as.character(catalog$reaction_id)), previous_union_ids
  )
  if (length(missing_from_union)) {
    stop(
      "Direct database cross-reference expansion produced reactions absent from the previous global union GEM: ",
      paste(utils::head(missing_from_union, 10L), collapse = ", "),
      ". Rebuild the original union with matching reaction annotations.",
      call. = FALSE
    )
  }
  union_match <- match(
    as.character(catalog$reaction_id),
    as.character(global_reaction_membership$reaction_id)
  )
  catalog$previous_union_is_core <- if (
    "is_core" %in% colnames(global_reaction_membership)
  ) {
    global_reaction_membership$is_core[union_match] %in% TRUE
  } else {
    as.character(catalog$reaction_id) %in%
      as.character(global_core_reactions$reaction_id)
  }
  catalog$previous_union_inclusion_stage <- if (
    "inclusion_stage" %in% colnames(global_reaction_membership)
  ) {
    as.character(global_reaction_membership$inclusion_stage[union_match])
  } else {
    NA_character_
  }
  catalog$score_target <- !catalog$previous_union_is_core
  catalog$target_role <- ifelse(
    catalog$previous_union_is_core,
    "previous_global_core_not_rescored",
    "direct_database_crossref_noncore"
  )
  catalog$lp_exclusion_reason <- ifelse(
    catalog$previous_union_is_core,
    "already_scored_in_original_layer2",
    NA_character_
  )
  target_relations <- catalog[catalog$score_target, , drop = FALSE]
  if (!nrow(target_relations)) {
    stop(
      "All directly linked KEGG, Reactome, or master-Rhea reactions were already scored as global cores in the original Layer 2 run.",
      call. = FALSE
    )
  }
  targets <- .rc_target_union_aggregate_targets(target_relations)
  rownames(selected) <- NULL
  rownames(catalog) <- NULL
  summary <- data.frame(
    n_selected_previous_core = nrow(selected),
    n_direct_crossref_relations = nrow(catalog),
    n_direct_crossref_reactions = length(unique(catalog$reaction_id)),
    n_previous_core_reactions_not_rescored = length(unique(
      catalog$reaction_id[catalog$previous_union_is_core]
    )),
    n_expanded_score_targets = nrow(targets),
    n_previous_union_reactions = length(previous_union_ids),
    gene_match = gene_match,
    expansion_policy =
      "direct_from_selected_core_via_kegg_reactome_master_rhea_only",
    scoring_policy =
      "direct_database_crossref_noncore_reactions_only",
    model_policy = "reuse_exact_previous_global_union_gem",
    stringsAsFactors = FALSE
  )
  list(
    selected_core_reactions = selected,
    expanded_reaction_catalog = catalog,
    expanded_scoring_targets = targets,
    previous_union_membership = global_reaction_membership,
    summary = summary,
    params = list(
      gene_match = gene_match,
      selected_core_reactions = unique(as.character(selected$reaction_id)),
      previous_core_reactions_not_rescored = unique(as.character(
        catalog$reaction_id[catalog$previous_union_is_core]
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
  if (!inherits(layer2, "regcompass_layer2_step") ||
      !identical(as.character(layer2$model_mode), "meta_module_gem")) {
    stop(
      "`layer2` must be the completed core LP stage with `model_mode = \"meta_module_gem\"`.",
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
  summary$medium_scenario <- trimws(as.character(summary$medium_scenario))
  summary$file <- as.character(summary$file)
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
      "Each medium scenario must resolve to one previous union GEM file: ",
      paste(ambiguous, collapse = ", "), call. = FALSE
    )
  }
  summary <- summary[!duplicated(summary$medium_scenario), , drop = FALSE]
  missing_files <- summary$file[!file.exists(summary$file)]
  if (length(missing_files)) {
    stop(
      "Previous union GEM cache files are unavailable: ",
      paste(utils::head(missing_files, 5L), collapse = ", "), call. = FALSE
    )
  }

  cache <- list()
  diagnostics <- list()
  fingerprints <- character(nrow(summary))
  for (i in seq_len(nrow(summary))) {
    scenario <- summary$medium_scenario[[i]]
    file <- summary$file[[i]]
    model <- readRDS(file)
    validated <- rc_validate_gem(model)
    fingerprints[[i]] <- .rc_full_gem_cache_fingerprint(model)
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
        module_id = "GLOBAL_UNION",
        reaction_id = reaction,
        target_direction = direction,
        medium_scenario = scenario,
        condition = "all",
        file = file,
        build_strategy = "reuse_exact_previous_global_union_gem"
      )
    }
  }
  if (!length(cache)) {
    stop("No direct database-linked non-core reaction direction can be scored.",
         call. = FALSE)
  }
  summary$source_model_fingerprint <- fingerprints
  summary$source_model_md5 <- unname(tools::md5sum(summary$file))
  summary$build_strategy <- "reuse_exact_previous_global_union_gem"
  summary$reused_without_rebuilding <- TRUE
  attr(cache, "summary") <- summary
  attr(cache, "direction_diagnostics") <- .rc_bind_frames_fill(diagnostics)
  cache
}

.rc_score_existing_union_cache <- function(
    layer1, gem, model_cache,
    condition_col, sample_col, celltype_col,
    omega = 0.95,
    solver = c("highs", "gurobi", "glpk"),
    time_limit = 60, flux_threshold = 1e-8,
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
    colnames(readRDS(model_cache[[row_id]]$file)$S)
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
    file = unique_files, unit_id = units,
    stringsAsFactors = FALSE
  )
  run_one <- function(task) {
    file <- as.character(task$file)
    unit_id <- as.character(task$unit_id)
    selected <- row_ids[model_files == file]
    model <- readRDS(file)
    answers <- lapply(selected, function(row_id) {
      entry <- model_cache[[row_id]]
      fit <- rc_compass_two_step_lp_directional(
        S = model$S, lb = model$lb, ub = model$ub,
        target_reaction = entry$reaction_id,
        penalties = penalties$penalty[colnames(model$S), unit_id],
        target_direction = entry$target_direction,
        omega = omega, solver = solver,
        time_limit = time_limit, flux_threshold = flux_threshold
      )
      list(
        row_id = row_id,
        unit_id = unit_id,
        penalty = fit$penalty,
        vmax = fit$vmax,
        feasible = isTRUE(fit$feasible),
        diagnostics = data.frame(
          row_id = row_id, unit_id = unit_id,
          sample_id = "global", module_id = "GLOBAL_UNION",
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
          source_union_model_file = file,
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
  model_diagnostics <- .rc_bind_frames_fill(lapply(seq_len(nrow(summary)), function(i) {
    model <- readRDS(summary$file[[i]])
    out <- model$closure_diagnostics %||% data.frame()
    if (nrow(out)) {
      out$medium_scenario <- summary$medium_scenario[[i]]
      out$source_union_model_file <- summary$file[[i]]
    }
    out
  }))
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
    md5 = unname(tools::md5sum(summary$file)),
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
    model_mode = "reused_global_union_gem",
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
      shared_gem_scope = "previous_global_union_by_medium",
      structural_model_reused_exactly = TRUE,
      parallel_task = "reused_union_model_by_metacell",
      flux_threshold = flux_threshold,
      solver = solver,
      time_limit = time_limit
    ),
    method = "microCOMPASS directional LP for direct KEGG/Reactome/master-Rhea-linked non-core reactions on exact previous global union GEMs"
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

#' Score directly database-linked non-core reactions in a previous union GEM
#'
#' Uses core reactions from a completed Layer 2 run as anchors. It directly
#' identifies reactions sharing KEGG, Reactome, or master-Rhea identifiers with
#' those anchors and scores only reactions that were not global core targets in
#' the original Layer 2 run. No subsystem or transitive expansion is performed.
#'
#' @param layer1 Output from [rc_regcompass_step_layer1()].
#' @param meta_modules Output from [rc_regcompass_step_meta_modules()].
#' @param layer2 Output from [rc_regcompass_step_layer2()] with
#'   `model_mode = "meta_module_gem"`.
#' @param gem The same GEM used for the original run.
#' @param outdir Output directory.
#' @param core_reaction_ids Previous core reaction IDs used as direct mapping
#'   anchors.
#' @param core_genes Genes used to resolve previous core anchors through GPRs.
#' @param gene_match Require a complete GPR group or allow any direct gene match.
#' @param layer2_args Optional `omega`, `target_direction`, `solver`,
#'   `time_limit`, and `flux_threshold` overrides.
#' @param parallel Whether to parallelize model-by-metacell tasks.
#' @param BPPARAM Optional BiocParallel parameter object.
#' @return A `regcompass_target_union_step` with selected core anchors, direct
#'   database relation rows, unique non-core LP targets, source-model provenance,
#'   and second-pass LP results.
#' @export
rc_regcompass_step_target_union <- function(
    layer1, meta_modules, layer2, gem, outdir,
    core_reaction_ids = NULL, core_genes = NULL,
    gene_match = c("complete_gpr", "any_direct"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL) {
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
  if (!is.list(layer2_args)) stop("`layer2_args` must be a list.", call. = FALSE)
  allowed <- c("omega", "target_direction", "solver", "time_limit", "flux_threshold")
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop("Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  global <- meta_modules$global_modules
  if (!is.list(global) ||
      !is.data.frame(global$global_core_reactions) ||
      !is.data.frame(global$global_reaction_membership)) {
    stop("The previous global meta-module union is unavailable.", call. = FALSE)
  }
  if (!setequal(
    as.character(layer2$source_core_reactions$reaction_id),
    as.character(global$global_core_reactions$reaction_id)
  )) {
    stop("Layer 2 was not generated from the supplied global core reactions.",
         call. = FALSE)
  }
  definition <- .rc_build_target_union_definition(
    gem = gem,
    global_core_reactions = global$global_core_reactions,
    global_reaction_membership = global$global_reaction_membership,
    core_reaction_ids = core_reaction_ids,
    core_genes = core_genes,
    gene_match = gene_match
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
  time_limit <- layer2_args$time_limit %||% 60
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
    time_limit = time_limit,
    flux_threshold = flux_threshold,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  scored$workflow_params <- workflow
  scored$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  scored$params$target_direction <- target_direction
  scored$params$target_scope <-
    "direct_kegg_reactome_master_rhea_noncore_only"
  scored$params$n_selected_previous_core <-
    nrow(definition$selected_core_reactions)
  scored$params$n_previous_core_reactions_not_rescored <-
    length(definition$params$previous_core_reactions_not_rescored)
  scored$params$n_expanded_score_targets <-
    length(definition$params$score_targets)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  .rc_mm_write_tsv_gz(
    definition$selected_core_reactions,
    file.path(outdir, "selected_previous_core_reactions.tsv.gz")
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
    definition$previous_union_membership,
    file.path(outdir, "reused_global_union_membership.tsv.gz")
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
  saveRDS(answer, file.path(outdir, "step_target_union.rds"))
  answer
}
