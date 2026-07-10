#' Build structural micro-GEMs once per target, direction, and medium scenario
#' @export
rc_build_microgem_cache <- function(gem, dirs, medium_scenarios, microgem_params = list()) {
  if (!is.data.frame(dirs) || !all(c("reaction_id", "target_direction") %in% colnames(dirs))) {
    stop("`dirs` must contain `reaction_id` and `target_direction`.", call. = FALSE)
  }
  if (is.null(medium_scenarios)) {
    medium_scenarios <- data.frame(
      medium_scenario_id = "base",
      exchange_reaction_id = character(),
      lb = numeric(),
      ub = numeric(),
      available = logical(),
      stringsAsFactors = FALSE
    )
  }
  if (!"medium_scenario_id" %in% colnames(medium_scenarios)) medium_scenarios$medium_scenario_id <- "custom"
  scenarios <- unique(as.character(medium_scenarios$medium_scenario_id))
  tasks <- expand.grid(dir_i = seq_len(nrow(dirs)), medium_scenario = scenarios, stringsAsFactors = FALSE)
  cache <- list()
  summary <- vector("list", nrow(tasks))
  strategy <- microgem_params$strategy %||% "target_khop"
  strict_closure <- isTRUE(microgem_params$strict_closure)
  closure_bad <- function(mg) {
    closure <- mg$closure_diagnostics
    n_support <- if (is.data.frame(closure) && "n_support_reactions" %in% colnames(closure) && nrow(closure)) closure$n_support_reactions[[1L]] else NA_real_
    target_feasible <- is.data.frame(closure) && "strict_target_feasible" %in% colnames(closure) && nrow(closure) && isTRUE(closure$strict_target_feasible[[1L]])
    !is.data.frame(closure) || !nrow(closure) ||
      !target_feasible ||
      (is.finite(n_support) && n_support == 0L)
  }
  for (k in seq_len(nrow(tasks))) {
    d <- dirs[tasks$dir_i[k], , drop = FALSE]
    sc <- tasks$medium_scenario[k]
    mt <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == sc, , drop = FALSE]
    key <- paste(d$reaction_id, d$target_direction, sc, sep = "::")
    params <- microgem_params
    params$strategy <- NULL
    params$module_col <- NULL
    if (identical(strategy, "auto")) params$strict_closure <- FALSE
    params <- params[names(params) %in% names(formals(rc_build_target_microgem))]
    mg <- do.call(rc_build_target_microgem,
                  c(list(gem = gem, target_reaction = d$reaction_id, medium_table = mt),
                    params))
    build_strategy <- "target_khop"
    fallback_from <- NA_character_
    target_status <- if (closure_bad(mg)) "structurally_infeasible" else "ok"
    if (identical(strategy, "auto") && closure_bad(mg)) {
      module_col <- microgem_params$module_col %||% "metabolic_module"
      module_params <- microgem_params
      module_params$strategy <- NULL
      module_params$module_col <- NULL
      module_params$strict_closure <- FALSE
      module_params <- module_params[names(module_params) %in% names(formals(rc_build_module_meso_gem))]
      mid <- .rc_module_id_for_reaction(gem, d$reaction_id, module_col = module_col)
      mg <- do.call(rc_build_module_meso_gem,
                    c(list(gem = gem, module_id = mid, medium_table = mt, module_col = module_col),
                      module_params))
      mg$target_reaction <- d$reaction_id
      mg$closure_diagnostics <- rc_check_module_gem_closure(mg, d$reaction_id)
      mg$build_params$fallback_from <- "target_khop"
      mg$build_params$strategy <- "auto"
      build_strategy <- "module_meso_gem"
      fallback_from <- "target_khop"
      target_status <- if (closure_bad(mg)) "structurally_infeasible" else "ok"
    }
    if (isTRUE(strict_closure) && closure_bad(mg)) {
      stop("Target micro-GEM failed strict closure validation: ", d$reaction_id, call. = FALSE)
    }
    mg$target_status <- target_status
    cache[[key]] <- mg
    summary[[k]] <- data.frame(
      cache_key = key,
      reaction_id = d$reaction_id,
      target_direction = d$target_direction,
      medium_scenario = sc,
      n_reactions = ncol(mg$S),
      n_metabolites = nrow(mg$S),
      k_hop = microgem_params$k_hop %||% 2,
      include_transport = microgem_params$include_transport %||% TRUE,
      include_exchange = microgem_params$include_exchange %||% TRUE,
      include_demand_sink = microgem_params$include_demand_sink %||% TRUE,
      include_cofactor_modules = microgem_params$include_cofactor_modules %||% TRUE,
      max_reactions = microgem_params$max_reactions %||% 500,
      build_strategy = build_strategy,
      fallback_from = fallback_from,
      target_status = target_status,
      model_version = (gem$model_info$model_version %||% gem$model_info$version %||% NA_character_),
      model_commit = (gem$model_info$source_commit %||% gem$model_info$commit %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }
  attr(cache, "summary") <- do.call(rbind, summary)
  cache
}
