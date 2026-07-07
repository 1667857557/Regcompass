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
  for (k in seq_len(nrow(tasks))) {
    d <- dirs[tasks$dir_i[k], , drop = FALSE]
    sc <- tasks$medium_scenario[k]
    mt <- medium_scenarios[as.character(medium_scenarios$medium_scenario_id) == sc, , drop = FALSE]
    key <- paste(d$reaction_id, d$target_direction, sc, sep = "::")
    mg <- do.call(rc_build_target_microgem,
                  c(list(gem = gem, target_reaction = d$reaction_id, medium_table = mt),
                    microgem_params))
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
      model_version = (gem$model_info$model_version %||% gem$model_info$version %||% NA_character_),
      model_commit = (gem$model_info$source_commit %||% gem$model_info$commit %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }
  attr(cache, "summary") <- do.call(rbind, summary)
  cache
}
