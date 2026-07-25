# Medium-specific union-GEM architecture.
#
# Stage 3 constructs biological meta-modules and a deduplicated merged catalogue.
# It does not run FASTCORE and does not create a GEM. Stage 5 constructs one
# medium-specific union GEM and performs the only FASTCORE completion.

.rc_build_condition_meta_modules <- function(grn_result, gem, outdir,
                                             layer1_args = list()) {
  if (!is.list(grn_result) ||
      !is.data.frame(grn_result$tf_peak_gene_significant)) {
    stop("`grn_result` is not a valid single-cell GRN result.", call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_cols <- grn_result$group_cols
  display_cols <- c("group_id", group_cols)
  module_cols <- unique(c(display_cols, "sample_id", "module_id"))
  metabolic_genes <- gem$metabolic_genes %||%
    rc_metabolic_gpr_genes(gem$gpr_table)
  sig <- grn_result$tf_peak_gene_significant
  sig$sample_id <- sig$group_id
  projection <- rc_project_metabolic_grn(
    sig,
    metabolic_genes = metabolic_genes,
    top_k = layer1_args$top_k_neighbors %||% 5L,
    min_shared_tfs = layer1_args$min_shared_tfs %||% 1L,
    min_tf_jaccard = layer1_args$min_tf_jaccard %||% 0,
    max_targets_per_tf = layer1_args$max_targets_per_tf %||% 200L,
    include_direct_metabolic_tf = TRUE
  )
  group_meta <- unique(
    grn_result$sample_status[, display_cols, drop = FALSE]
  )
  group_meta$analysis_unit_id <- group_meta$group_id
  projection$nodes <- .rc_remap_projection_metadata(
    projection$nodes, group_meta, "analysis_unit_id", display_cols
  )
  projection$edges <- .rc_remap_projection_metadata(
    projection$edges, group_meta, "analysis_unit_id", display_cols
  )
  core <- rc_map_meta_module_core_reactions(
    projection$nodes, gem$gpr_table
  )
  if (nrow(core)) {
    core <- merge(
      core,
      unique(projection$nodes[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    core <- core[, c(
      display_cols, setdiff(colnames(core), display_cols)
    ), drop = FALSE]
  }
  expanded <- rc_expand_meta_module_reactions(
    gem,
    core,
    subsystem_table = layer1_args$subsystem_table %||% NULL,
    expansion_mode = layer1_args$expansion_mode %||% "ordered_once"
  )
  if (nrow(expanded$reaction_membership)) {
    expanded$reaction_membership <- merge(
      expanded$reaction_membership,
      unique(core[, module_cols, drop = FALSE]),
      by = c("sample_id", "module_id"),
      all.x = TRUE,
      sort = FALSE
    )
    expanded$reaction_membership <- expanded$reaction_membership[, c(
      display_cols,
      setdiff(colnames(expanded$reaction_membership), display_cols)
    ), drop = FALSE]
  }
  out <- c(grn_result, list(
    metabolic_gene_nodes = projection$nodes,
    metabolic_gene_edges = projection$edges,
    core_gene_reaction = core,
    reaction_membership = expanded$reaction_membership,
    biological_reaction_membership = expanded$reaction_membership,
    meta_module_summary = expanded$summary,
    crossref_maps = expanded$crossref_maps,
    analysis_group_unit = "condition_x_celltype_single_cell_grn",
    feasibility_completion = "none_at_meta_module_stage"
  ))
  .rc_mm_write_tsv_gz(
    projection$nodes, file.path(outdir, "metabolic_gene_nodes.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    projection$edges, file.path(outdir, "metabolic_gene_edges.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    core, file.path(outdir, "core_gene_reaction.tsv.gz")
  )
  .rc_mm_write_tsv_gz(
    expanded$reaction_membership,
    file.path(outdir, "meta_module_reactions.tsv.gz")
  )
  saveRDS(out, file.path(outdir, "condition_meta_modules.rds"))
  out
}

.rc_merge_stratum_meta_modules <- function(artifacts) {
  names_to_merge <- c(
    "sample_status", "tf_peak_gene_all", "tf_peak_gene_significant",
    "metabolic_gene_nodes", "metabolic_gene_edges", "core_gene_reaction",
    "reaction_membership", "meta_module_summary"
  )
  out <- lapply(names_to_merge, function(name) {
    .rc_bind_frames_fill(lapply(
      artifacts,
      function(artifact) artifact$grn_meta_modules[[name]]
    ))
  })
  names(out) <- names_to_merge

  core <- out$core_gene_reaction
  if ("is_core" %in% colnames(core)) {
    core <- core[core$is_core %in% TRUE, , drop = FALSE]
  }
  core_ids <- unique(as.character(core$reaction_id))
  core_ids <- core_ids[!is.na(core_ids) & nzchar(core_ids)]

  biological <- out$reaction_membership
  biological_ids <- unique(as.character(biological$reaction_id))
  biological_ids <- biological_ids[
    !is.na(biological_ids) & nzchar(biological_ids)
  ]
  if (!length(core_ids) || !length(biological_ids)) {
    stop("No merged biological meta-module reactions were produced.",
         call. = FALSE)
  }

  out$biological_reaction_membership <- biological
  out$merged_core_reactions <- data.frame(
    sample_id = "merged",
    module_id = "MERGED_META_MODULES",
    reaction_id = core_ids,
    is_core = TRUE,
    stringsAsFactors = FALSE
  )
  out$merged_reaction_membership <- data.frame(
    sample_id = "merged",
    module_id = "MERGED_META_MODULES",
    reaction_id = biological_ids,
    is_core = biological_ids %in% core_ids,
    inclusion_stage = ifelse(
      biological_ids %in% core_ids,
      "merged_meta_module_core",
      "merged_meta_module_biological_member"
    ),
    stringsAsFactors = FALSE
  )
  out$schema_version <- "regcompass_merged_meta_modules_v1"
  out$source_group_ids <- unique(vapply(
    artifacts,
    function(artifact) as.character(artifact$group_id),
    character(1)
  ))
  out$merge_source <- "deduplicated_biological_meta_module_reactions"
  out$is_gem <- FALSE
  out$fastcore_applied <- FALSE
  out
}

.rc_build_medium_specific_union_gem_cache <- function(
    gem, reaction_membership, core_reactions,
    target_reactions = NULL, medium_scenarios = NULL,
    cache_dir = tempfile("RegCompassR_medium_union_gem_"),
    target_direction = c("both", "forward", "reverse"),
    solver = "highs", time_limit = 300,
    fastcore_epsilon = 1e-4,
    max_support_reactions = 2000,
    strict = TRUE) {
  target_direction <- match.arg(target_direction)
  required <- c("sample_id", "module_id", "reaction_id")
  if (!is.data.frame(reaction_membership) ||
      !all(required %in% colnames(reaction_membership))) {
    stop("`reaction_membership` must contain sample_id, module_id and reaction_id.",
         call. = FALSE)
  }
  if (!is.data.frame(core_reactions) ||
      !all(required %in% colnames(core_reactions))) {
    stop("`core_reactions` must contain sample_id, module_id and reaction_id.",
         call. = FALSE)
  }
  if ("is_core" %in% colnames(core_reactions)) {
    core_reactions <- core_reactions[
      core_reactions$is_core %in% TRUE, , drop = FALSE
    ]
  }
  if (!is.null(target_reactions)) {
    core_reactions <- core_reactions[
      as.character(core_reactions$reaction_id) %in%
        as.character(target_reactions),
      , drop = FALSE
    ]
  }
  if (!nrow(core_reactions)) {
    stop("No merged core reactions remain for union-GEM scoring.",
         call. = FALSE)
  }
  reaction_membership$sample_id <- "global"
  reaction_membership$module_id <- "MEDIUM_UNION_GEM"
  core_reactions$sample_id <- "global"
  core_reactions$module_id <- "MEDIUM_UNION_GEM"
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  scenarios <- unique(as.character(
    medium_scenarios$medium_scenario_id
  ))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  safe <- function(value) {
    paste(sprintf("%02x", as.integer(charToRaw(enc2utf8(value)))),
          collapse = "")
  }
  cache <- list()
  summaries <- list()
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
    model <- rc_build_meta_module_gem(
      gem = gem,
      reaction_membership = reaction_membership,
      core_reactions = core_reactions,
      sample_id = "global",
      module_id = "MEDIUM_UNION_GEM",
      medium_table = medium,
      condition = NULL,
      target_direction = target_direction,
      solver = solver,
      time_limit = time_limit,
      fastcore_epsilon = fastcore_epsilon,
      max_support_reactions = max_support_reactions,
      strict = strict
    )
    model$shared_across_units <- TRUE
    model$is_union_gem <- TRUE
    model$union_gem_scope <-
      "one_medium_shared_across_conditions_and_metacells"
    model$union_gem_medium_scenario <- scenario
    model$build_params$strategy <- "medium_specific_union_gem"
    model$build_params$completion_stage <-
      "single_global_fastcore_after_meta_module_merge"
    file <- file.path(
      cache_dir,
      paste0("union_gem__medium_", safe(scenario), ".rds")
    )
    saveRDS(model, file)
    summaries[[scenario]] <- data.frame(
      medium_scenario = scenario,
      file = file,
      n_reactions = ncol(model$S),
      n_metabolites = nrow(model$S),
      n_biological_reactions =
        model$build_params$n_biological_reactions,
      n_fastcore_support_reactions =
        model$build_params$n_fastcore_support_reactions,
      n_merged_biological_reactions =
        model$build_params$n_biological_reactions,
      n_global_fastcore_support_reactions =
        model$build_params$n_fastcore_support_reactions,
      target_status = model$target_status,
      build_strategy = "medium_specific_union_gem",
      stringsAsFactors = FALSE
    )
    if (!nrow(model$target_directions)) next
    for (i in seq_len(nrow(model$target_directions))) {
      reaction <- as.character(
        model$target_directions$reaction_id[[i]]
      )
      direction <- as.character(
        model$target_directions$target_direction[[i]]
      )
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
        build_strategy = "medium_specific_union_gem"
      )
    }
  }
  attr(cache, "summary") <- .rc_bind_frames_fill(summaries)
  cache
}

# Internal bridge for the existing microCOMPASS engine. This symbol is not
# exported; the returned artefacts and public documentation use union-GEM terms.
.rc_build_global_meta_module_gem_cache <- function(...) {
  .rc_build_medium_specific_union_gem_cache(...)
}

#' Construct core reactions and biological meta-modules from GRNs
#' @export
rc_regcompass_step_meta_modules <- function(
    grn, metacells, gem, outdir,
    layer1_args = list(),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("meta_modules", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  if (!identical(.rc_workflow_signature(grn),
                 .rc_workflow_signature(metacells))) {
    stop("GRN and metacell stages use different metadata or assay settings.",
         call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  validated_gem <- rc_validate_gem(gem)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  group_coverage <- .rc_validate_grn_metacell_group_coverage(
    grn_result = grn$grn_result,
    metacell_meta = metacells$pooled$metacell_meta,
    condition_col = metacells$params$condition_col,
    celltype_col = metacells$params$celltype_col
  )
  .rc_write_tsv_gz(
    group_coverage,
    file.path(outdir, "grn_metacell_group_coverage.tsv.gz")
  )
  condition_modules <- .rc_build_condition_meta_modules(
    grn$grn_result, gem, outdir, layer1_args
  )
  condition_modules$grn_metacell_group_coverage <- group_coverage
  if (!is.data.frame(condition_modules$reaction_membership) ||
      !nrow(condition_modules$reaction_membership)) {
    stop("Meta-module construction produced no reaction membership.",
         call. = FALSE)
  }
  missing <- setdiff(
    unique(as.character(condition_modules$reaction_membership$reaction_id)),
    colnames(validated_gem$S)
  )
  if (length(missing)) {
    stop("Meta-module reactions absent from the GEM: ",
         paste(utils::head(missing, 10L), collapse = ", "),
         call. = FALSE)
  }
  merged_modules <- .rc_merge_stratum_meta_modules(list(list(
    group_id = "condition_pooled",
    grn_meta_modules = condition_modules
  )))
  if (!is.data.frame(merged_modules$merged_core_reactions) ||
      !nrow(merged_modules$merged_core_reactions)) {
    stop("No complete-GPR merged core reactions remain after module merging.",
         call. = FALSE)
  }
  answer <- list(
    condition_modules = condition_modules,
    merged_modules = merged_modules,
    group_coverage = group_coverage,
    workflow_params = metacells$params,
    grn_params = grn$params,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      layer1_args = layer1_args,
      local_fastcore = FALSE,
      merge_creates_gem = FALSE
    )
  )
  class(answer) <- c("regcompass_meta_module_step", "list")
  answer <- .rc_step_monitor_finish(answer, monitor)
  saveRDS(
    condition_modules,
    file.path(outdir, "condition_meta_modules.rds")
  )
  saveRDS(
    merged_modules,
    file.path(outdir, "merged_meta_modules.rds")
  )
  saveRDS(answer, file.path(outdir, "step_meta_modules.rds"))
  answer
}

#' Run directional COMPASS-like LP scoring
#' @export
rc_regcompass_step_layer2 <- function(
    layer1, meta_modules, gem, medium_scenarios, outdir,
    model_mode = c("meta_module_gem", "full_gem"),
    layer2_args = list(), parallel = TRUE, BPPARAM = NULL,
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("layer2", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  model_mode <- match.arg(model_mode)
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  if (!is.list(layer2_args)) {
    stop("`layer2_args` must be a list.", call. = FALSE)
  }
  params <- meta_modules$workflow_params
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  medium_scenarios <- .rc_validate_shared_medium(medium_scenarios)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  layer2_args$model_params <- layer2_args$model_params %||% list()
  layer2_args$model_params$cache_dir <- file.path(
    outdir, "model_cache", model_mode
  )
  reserved <- intersect(names(layer2_args), c(
    "layer1", "gem", "mode", "unit", "reaction_membership",
    "core_reactions", "target_reactions", "medium_scenarios",
    "sample_col", "condition_col", "celltype_col", "BPPARAM",
    "parallel", "penalty_weights"
  ))
  if (length(reserved)) {
    stop("`layer2_args` cannot override workflow fields: ",
         paste(reserved, collapse = ", "), call. = FALSE)
  }
  solver <- match.arg(
    as.character(layer2_args$solver %||% "highs"),
    c("highs", "gurobi", "glpk")
  )
  .rc_require_lp_solver(solver)
  catalogue <- meta_modules$merged_modules
  targets <- unique(as.character(
    catalogue$merged_core_reactions$reaction_id
  ))
  missing_expression <- setdiff(
    targets, rownames(layer1$reaction_expression)
  )
  if (length(missing_expression)) {
    stop("Merged core reactions are absent from Layer 1 expression: ",
         paste(utils::head(missing_expression, 10L), collapse = ", "),
         call. = FALSE)
  }
  defaults <- list(
    layer1 = layer1,
    gem = gem,
    target_reactions = targets,
    medium_scenarios = medium_scenarios,
    mode = model_mode,
    reaction_membership = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_reaction_membership
    } else {
      NULL
    },
    core_reactions = if (identical(model_mode, "meta_module_gem")) {
      catalogue$merged_core_reactions
    } else {
      NULL
    },
    unit = "metacell",
    sample_col = params$sample_col,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col,
    parallel = parallel,
    BPPARAM = BPPARAM
  )
  defaults[names(layer2_args)] <- NULL
  answer <- withCallingHandlers(
    do.call(rc_run_microcompass, c(defaults, layer2_args)),
    warning = function(w) {
      if (grepl(
        "Metacell-level scores are descriptive pseudo-observations",
        conditionMessage(w), fixed = TRUE
      )) invokeRestart("muffleWarning")
    }
  )
  answer$workflow_params <- params
  answer$gem_fingerprint <- .rc_stage_gem_fingerprint(gem)
  answer$source_core_reactions <- catalogue$merged_core_reactions
  answer$source_merged_reaction_membership <-
    catalogue$merged_reaction_membership
  answer$union_gem_policy <- if (identical(model_mode, "meta_module_gem")) {
    "one medium-specific union GEM; single global FASTCORE completion"
  } else {
    "shared full GEM; no union-GEM reconstruction"
  }
  class(answer) <- c("regcompass_layer2_step", "list")
  .rc_validate_layer2_stage(
    answer, layer1 = layer1, workflow_params = params, gem = gem,
    argument = "layer2"
  )
  answer <- .rc_step_monitor_finish(answer, monitor)
  rc_export_microcompass(answer, outdir)
  saveRDS(answer, file.path(outdir, "step_layer2.rds"))
  answer
}

#' Assemble final RegCompass results
#' @export
rc_regcompass_step_results <- function(
    grn, metacells, meta_modules, layer1, layer2, gem, outdir,
    species = c("auto", "human", "mouse"),
    progress = getOption("RegCompassR.progress", TRUE)) {
  monitor <- .rc_step_monitor_start("results", outdir, progress)
  on.exit(.rc_step_monitor_fail(monitor), add = TRUE)
  .rc_require_stage_class(
    grn, "regcompass_grn_step", "grn", "rc_regcompass_step_grn"
  )
  .rc_require_stage_class(
    metacells, "regcompass_metacell_step", "metacells",
    "rc_regcompass_step_metacells"
  )
  .rc_require_stage_class(
    meta_modules, "regcompass_meta_module_step", "meta_modules",
    "rc_regcompass_step_meta_modules"
  )
  params <- metacells$params
  if (!identical(params, meta_modules$workflow_params) ||
      !identical(.rc_workflow_signature(grn),
                 .rc_workflow_signature(metacells))) {
    stop("Upstream stages use different workflow parameters.",
         call. = FALSE)
  }
  .rc_require_stage_gem(grn, gem, "grn")
  .rc_require_stage_gem(meta_modules, gem, "meta_modules")
  .rc_validate_layer1_stage(
    layer1, workflow_params = params, gem = gem, argument = "layer1"
  )
  .rc_validate_layer2_stage(
    layer2, layer1 = layer1, workflow_params = params, gem = gem,
    argument = "layer2"
  )
  species <- .rc_infer_gem_species(gem, species)
  comparison <- .rc_condition_penalty_comparison(
    layer2,
    condition_col = params$condition_col,
    celltype_col = params$celltype_col
  )
  conditions <- unique(as.character(
    metacells$pooled$metacell_meta[[params$condition_col]]
  ))
  condition_fields <- intersect(c(
    "metabolic_gene_nodes", "metabolic_gene_edges",
    "core_gene_reaction", "biological_reaction_membership",
    "reaction_membership", "meta_module_summary",
    "analysis_group_unit", "grn_metacell_group_coverage",
    "feasibility_completion"
  ), names(meta_modules$condition_modules))
  condition_modules <-
    meta_modules$condition_modules[condition_fields]
  result <- list(
    schema_version = "regcompass_grn_first_v3",
    version = "1.8.3",
    species = species,
    model_mode = layer2$model_mode,
    analysis_mode = comparison$analysis_mode,
    grn = grn$grn_result,
    metacells = metacells$pooled,
    layer1 = layer1,
    condition_grn_meta_modules = condition_modules,
    merged_grn_meta_modules = meta_modules$merged_modules,
    grn_meta_modules = meta_modules$merged_modules,
    grn_metacell_group_coverage = meta_modules$group_coverage,
    microcompass = layer2,
    reaction_ranking = comparison$ranking,
    condition_summary = comparison$summary,
    condition_contrast = comparison$contrast,
    inference_policy = comparison$inference_policy,
    gem_fingerprint = .rc_stage_gem_fingerprint(gem),
    params = list(
      n_conditions = length(conditions),
      workflow_order = c(
        "single_cell_grn", "condition_metacells", "meta_modules",
        "layer1", "medium_specific_union_gem_layer2"
      ),
      pando_grouping = c(params$condition_col, params$celltype_col),
      pando_peak_cor =
        grn$grn_result$normalization_policy$pando_peak_cor,
      metacell_grouping = params$condition_col,
      metacell_celltype_assignment =
        "supercell_label_guided_then_dominant_membership_audit",
      metacell_gamma = params$metacell_args$gamma,
      sample_weighting = "none",
      meta_module_expansion =
        "core_subsystem_plus_kegg_reactome_master_rhea_only",
      meta_module_merge =
        "reaction_id_deduplication_only_not_a_gem",
      feasibility_completion = if (
        identical(layer2$model_mode, "meta_module_gem")
      ) {
        "single_global_fastcore_on_each_medium_specific_union_gem"
      } else {
        "not_applicable_full_gem"
      },
      feasibility_completion_stages = if (
        identical(layer2$model_mode, "meta_module_gem")
      ) {
        "layer2_medium_specific_only"
      } else {
        "none"
      },
      union_gem_definition =
        "medium_constrained_merged_meta_modules_plus_global_fastcore_support",
      pando_normalization_policy =
        grn$grn_result$normalization_policy,
      penalty_formula = "1/(1+log2(1+E_multiome))",
      execution_mode = "stepwise"
    )
  )
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
  result <- .rc_step_monitor_finish(result, monitor)
  saveRDS(comparison, file.path(outdir, "step_comparison.rds"))
  saveRDS(result, file.path(outdir, "regcompass_result.rds"))
  result
}

#' Score directly database-linked non-core reactions in previous union GEMs
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
  allowed <- c(
    "omega", "target_direction", "solver", "time_limit",
    "flux_threshold"
  )
  unknown <- setdiff(names(layer2_args), allowed)
  if (length(unknown)) {
    stop("Unsupported `layer2_args`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  catalogue <- meta_modules$merged_modules
  if (!is.list(catalogue) ||
      !is.data.frame(catalogue$merged_core_reactions) ||
      !is.data.frame(catalogue$merged_reaction_membership)) {
    stop("The merged meta-module catalogue is unavailable.",
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
    global_core_reactions = catalogue$merged_core_reactions,
    global_reaction_membership = catalogue$merged_reaction_membership,
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
  scored$params$n_cached_union_unavailable_reactions <-
    length(definition$params$cached_union_unavailable_reactions)
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
    file.path(outdir, "reused_union_gem_membership.tsv.gz")
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
