# GRN-to-reaction biological meta-module definition.

#' Map GRN metabolic gene modules to core Human-GEM reactions
rc_map_meta_module_core_reactions <- function(gene_nodes, gpr_table) {
  if (!is.data.frame(gene_nodes) ||
      !all(c("sample_id", "gene", "module_id") %in% colnames(gene_nodes))) {
    stop("`gene_nodes` must contain sample_id, gene and module_id.",
         call. = FALSE)
  }
  if (!is.data.frame(gpr_table) ||
      !all(c("reaction_id", "and_group_id", "gene") %in% colnames(gpr_table))) {
    stop(
      "`gpr_table` must contain reaction_id, and_group_id and gene; treating every gene as an isoenzyme would misclassify required subunits.",
      call. = FALSE
    )
  }

  nodes <- unique(gene_nodes[, c("sample_id", "module_id", "gene"), drop = FALSE])
  nodes$sample_id <- as.character(nodes$sample_id)
  nodes$module_id <- as.character(nodes$module_id)
  nodes$gene <- toupper(trimws(as.character(nodes$gene)))
  nodes <- nodes[!is.na(nodes$gene) & nzchar(nodes$gene), , drop = FALSE]

  gpr <- gpr_table
  gpr$reaction_id <- trimws(as.character(gpr$reaction_id))
  gpr$gene <- toupper(trimws(as.character(gpr$gene)))
  gpr <- gpr[
    !is.na(gpr$reaction_id) & nzchar(gpr$reaction_id) &
      !is.na(gpr$gene) & nzchar(gpr$gene),
    , drop = FALSE
  ]
  gpr$and_group_id <- as.character(gpr$and_group_id)
  gpr <- unique(gpr[, c("reaction_id", "and_group_id", "gene"), drop = FALSE])

  empty <- data.frame(
    sample_id = character(),
    module_id = character(),
    gene = character(),
    reaction_id = character(),
    and_group_id = character(),
    required_genes = character(),
    matched_genes = character(),
    missing_genes = character(),
    group_complete = logical(),
    is_core = logical(),
    is_partial_candidate = logical(),
    inclusion_stage = character(),
    stringsAsFactors = FALSE
  )
  if (!nrow(nodes) || !nrow(gpr)) return(empty)

  module_rows <- split(
    seq_len(nrow(nodes)),
    paste(nodes$sample_id, nodes$module_id, sep = "\001")
  )
  output <- list()
  output_index <- 0L

  for (module_key in names(module_rows)) {
    module_index <- module_rows[[module_key]]
    module_genes <- unique(nodes$gene[module_index])
    sample_id <- nodes$sample_id[module_index[[1L]]]
    module_id <- nodes$module_id[module_index[[1L]]]
    candidate_reactions <- unique(gpr$reaction_id[gpr$gene %in% module_genes])

    for (reaction_id in candidate_reactions) {
      reaction_gpr <- gpr[gpr$reaction_id == reaction_id, , drop = FALSE]
      group_rows <- split(seq_len(nrow(reaction_gpr)), reaction_gpr$and_group_id)
      group_complete <- vapply(group_rows, function(rows) {
        all(unique(reaction_gpr$gene[rows]) %in% module_genes)
      }, logical(1))
      reaction_core <- any(group_complete)

      for (group_id in names(group_rows)) {
        rows <- group_rows[[group_id]]
        required <- sort(unique(reaction_gpr$gene[rows]))
        matched <- intersect(required, module_genes)
        if (!length(matched)) next
        missing <- setdiff(required, module_genes)
        complete <- !length(missing)

        for (gene in matched) {
          output_index <- output_index + 1L
          output[[output_index]] <- data.frame(
            sample_id = sample_id,
            module_id = module_id,
            gene = gene,
            reaction_id = reaction_id,
            and_group_id = group_id,
            required_genes = paste(required, collapse = ";"),
            matched_genes = paste(sort(matched), collapse = ";"),
            missing_genes = paste(sort(missing), collapse = ";"),
            group_complete = complete,
            is_core = reaction_core,
            is_partial_candidate = !reaction_core,
            inclusion_stage = if (reaction_core) {
              "core_complete_gpr"
            } else {
              "partial_gpr_candidate"
            },
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (!length(output)) return(empty)
  answer <- unique(do.call(rbind, output))
  rownames(answer) <- NULL
  answer
}

.rc_clean_meta_module_map <- function(x, id_col) {
  if (!is.data.frame(x) ||
      !all(c("reaction_id", id_col) %in% colnames(x))) {
    return(data.frame())
  }
  reaction <- trimws(as.character(x$reaction_id))
  identifier <- trimws(as.character(x[[id_col]]))
  keep <- !is.na(reaction) & nzchar(reaction) &
    !is.na(identifier) & nzchar(identifier)
  out <- x[keep, c("reaction_id", id_col), drop = FALSE]
  out$reaction_id <- reaction[keep]
  out[[id_col]] <- identifier[keep]
  unique(out)
}

#' Expand core reactions into GRN-defined biological reaction meta-modules
#'
#' Expansion is restricted to reactions in core-reaction subsystems and reactions
#' sharing KEGG, Reactome or master-Rhea identifiers. No stoichiometric-neighbour
#' expansion is performed. Flux-feasibility support is added later by FASTCORE.
.rc_expand_meta_module_reactions_core <- function(
    gem, core_reactions, subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L) {
  expansion_mode <- match.arg(expansion_mode)
  if (!is.numeric(max_iterations) || length(max_iterations) != 1L ||
      !is.finite(max_iterations) || max_iterations < 1 ||
      abs(max_iterations - round(max_iterations)) > sqrt(.Machine$double.eps)) {
    stop("`max_iterations` must be one positive integer.", call. = FALSE)
  }
  required <- c("sample_id", "module_id", "gene", "reaction_id")
  if (!is.data.frame(core_reactions) ||
      !all(required %in% colnames(core_reactions))) {
    stop(
      paste(
        "`core_reactions` must contain sample_id, module_id,",
        "gene and reaction_id."
      ),
      call. = FALSE
    )
  }

  maps <- rc_reaction_crossref_maps(gem, subsystem_table = subsystem_table)
  maps$subsystem <- .rc_clean_meta_module_map(maps$subsystem, "subsystem_id")
  maps$kegg <- .rc_clean_meta_module_map(maps$kegg, "kegg_id")
  maps$reactome <- .rc_clean_meta_module_map(maps$reactome, "reactome_id")
  maps$rhea_master <- .rc_clean_meta_module_map(
    maps$rhea_master, "rhea_master_id"
  )
  if (!nrow(maps$subsystem)) {
    stop("No usable reaction-to-subsystem annotations were found.", call. = FALSE)
  }

  validated <- rc_validate_gem(gem)
  valid_reactions <- colnames(validated$S)

  if ("is_core" %in% colnames(core_reactions)) {
    core_reactions <- core_reactions[
      core_reactions$is_core %in% TRUE,
      , drop = FALSE
    ]
  }
  core_reactions <- core_reactions[
    core_reactions$reaction_id %in% valid_reactions,
    , drop = FALSE
  ]
  if (!nrow(core_reactions)) {
    stop("No GRN genes mapped to valid GEM reactions.", call. = FALSE)
  }

  groups <- split(
    seq_len(nrow(core_reactions)),
    paste(core_reactions$sample_id, core_reactions$module_id, sep = "\001")
  )
  membership_rows <- list()
  summary_rows <- list()

  for (key in names(groups)) {
    index <- groups[[key]]
    core <- unique(as.character(core_reactions$reaction_id[index]))
    sample_id <- as.character(core_reactions$sample_id[index[[1L]]])
    module_id <- as.character(core_reactions$module_id[index[[1L]]])
    included <- core
    reasons <- stats::setNames(rep("core_grn_gene", length(core)), core)
    source_ids <- stats::setNames(rep(NA_character_, length(core)), core)

    add_reactions <- function(reactions, reason, source_map = NULL) {
      reactions <- intersect(.rc_mm_trim_unique(reactions), valid_reactions)
      new <- setdiff(reactions, included)
      if (length(new)) {
        included <<- c(included, new)
        reasons[new] <<- reason
        if (!is.null(source_map)) source_ids[new] <<- source_map[new]
      }
      invisible(new)
    }
    reactions_in_subsystems <- function(subsystems) {
      unique(maps$subsystem$reaction_id[
        maps$subsystem$subsystem_id %in% subsystems
      ])
    }
    subsystems_for_reactions <- function(reactions) {
      unique(maps$subsystem$subsystem_id[
        maps$subsystem$reaction_id %in% reactions
      ])
    }
    source_for_subsystem_members <- function(subsystems) {
      reactions <- reactions_in_subsystems(subsystems)
      stats::setNames(vapply(reactions, function(reaction) {
        ids <- unique(maps$subsystem$subsystem_id[
          maps$subsystem$reaction_id == reaction &
            maps$subsystem$subsystem_id %in% subsystems
        ])
        paste(paste0("SUBSYSTEM:", ids), collapse = ";")
      }, character(1)), reactions)
    }
    prefixed_crossrefs <- function(prefix, values) {
      values <- .rc_mm_trim_unique(values)
      if (!length(values)) character() else paste0(prefix, values)
    }

    expand_annotations_once <- function() {
      before <- length(included)
      core_subsystems <- subsystems_for_reactions(core)
      add_reactions(
        reactions_in_subsystems(core_subsystems),
        "same_core_subsystem",
        source_for_subsystem_members(core_subsystems)
      )

      current <- included
      kegg_ids <- unique(maps$kegg$kegg_id[
        maps$kegg$reaction_id %in% current
      ])
      reactome_ids <- unique(maps$reactome$reactome_id[
        maps$reactome$reaction_id %in% current
      ])
      database_reactions <- unique(c(
        maps$kegg$reaction_id[maps$kegg$kegg_id %in% kegg_ids],
        maps$reactome$reaction_id[
          maps$reactome$reactome_id %in% reactome_ids
        ]
      ))
      database_source <- stats::setNames(
        vapply(database_reactions, function(reaction) {
          reaction_kegg <- intersect(
            unique(maps$kegg$kegg_id[
              maps$kegg$reaction_id == reaction
            ]),
            kegg_ids
          )
          reaction_reactome <- intersect(
            unique(maps$reactome$reactome_id[
              maps$reactome$reaction_id == reaction
            ]),
            reactome_ids
          )
          paste(
            c(
              prefixed_crossrefs("KEGG:", reaction_kegg),
              prefixed_crossrefs("REACTOME:", reaction_reactome)
            ),
            collapse = ";"
          )
        }, character(1)),
        database_reactions
      )
      add_reactions(
        database_reactions,
        "shared_kegg_or_reactome_reaction",
        database_source
      )

      current <- included
      master_ids <- unique(maps$rhea_master$rhea_master_id[
        maps$rhea_master$reaction_id %in% current
      ])
      rhea_reactions <- unique(maps$rhea_master$reaction_id[
        maps$rhea_master$rhea_master_id %in% master_ids
      ])
      rhea_source <- stats::setNames(
        vapply(rhea_reactions, function(reaction) {
          ids <- intersect(
            unique(maps$rhea_master$rhea_master_id[
              maps$rhea_master$reaction_id == reaction
            ]),
            master_ids
          )
          paste(paste0("RHEA_MASTER:", ids), collapse = ";")
        }, character(1)),
        rhea_reactions
      )
      add_reactions(
        rhea_reactions,
        "shared_master_rhea_reaction",
        rhea_source
      )
      length(included) > before
    }

    changed <- expand_annotations_once()
    iteration <- 1L
    if (identical(expansion_mode, "fixed_point")) {
      while (isTRUE(changed) && iteration < as.integer(max_iterations)) {
        iteration <- iteration + 1L
        changed <- expand_annotations_once()
      }
    }

    membership <- data.frame(
      sample_id = sample_id,
      module_id = module_id,
      reaction_id = included,
      is_core = included %in% core,
      inclusion_stage = unname(reasons[included]),
      source_annotation = unname(source_ids[included]),
      expansion_mode = expansion_mode,
      stringsAsFactors = FALSE
    )
    membership_rows[[key]] <- membership
    summary_rows[[key]] <- data.frame(
      sample_id = sample_id,
      module_id = module_id,
      n_core_genes = length(unique(core_reactions$gene[index])),
      n_core_reactions = length(core),
      n_reactions = length(included),
      n_subsystem_added = sum(
        membership$inclusion_stage == "same_core_subsystem"
      ),
      n_database_added = sum(
        membership$inclusion_stage ==
          "shared_kegg_or_reactome_reaction"
      ),
      n_rhea_added = sum(
        membership$inclusion_stage ==
          "shared_master_rhea_reaction"
      ),
      iterations = iteration,
      stringsAsFactors = FALSE
    )
  }

  list(
    reaction_membership = do.call(rbind, membership_rows),
    summary = do.call(rbind, summary_rows),
    crossref_maps = maps
  )
}

rc_expand_meta_module_reactions <- function(
    gem, core_reactions, subsystem_table = NULL,
    expansion_mode = c("ordered_once", "fixed_point"),
    max_iterations = 10L) {
  hard_core_reactions <- .rc_hard_core_rows(core_reactions)
  answer <- .rc_expand_meta_module_reactions_core(
    gem = gem,
    core_reactions = hard_core_reactions,
    subsystem_table = subsystem_table,
    expansion_mode = expansion_mode,
    max_iterations = max_iterations
  )
  if (!"is_core" %in% colnames(core_reactions)) return(answer)

  key <- function(data) paste(
    as.character(data$sample_id),
    as.character(data$module_id),
    as.character(data$reaction_id),
    sep = "\001"
  )
  candidate_keys <- unique(key(core_reactions))
  hard_rows <- .rc_hard_core_rows(core_reactions)
  hard_keys <- unique(key(hard_rows))
  membership_keys <- key(answer$reaction_membership)

  answer$reaction_membership$is_core <- membership_keys %in% hard_keys
  partial_anchor <- membership_keys %in% setdiff(candidate_keys, hard_keys)
  if (any(partial_anchor)) {
    answer$reaction_membership <- answer$reaction_membership[
      !partial_anchor,
      , drop = FALSE
    ]
  }

  if (is.data.frame(answer$summary) && nrow(answer$summary)) {
    for (i in seq_len(nrow(answer$summary))) {
      sample_id <- as.character(answer$summary$sample_id[[i]])
      module_id <- as.character(answer$summary$module_id[[i]])
      selected <- as.character(core_reactions$sample_id) == sample_id &
        as.character(core_reactions$module_id) == module_id &
        core_reactions$is_core %in% TRUE
      gene_selected <- selected
      if ("group_complete" %in% colnames(core_reactions)) {
        gene_selected <- gene_selected &
          core_reactions$group_complete %in% TRUE
      }

      membership_selected <-
        as.character(answer$reaction_membership$sample_id) == sample_id &
        as.character(answer$reaction_membership$module_id) == module_id
      membership <- answer$reaction_membership[
        membership_selected,
        , drop = FALSE
      ]
      count_stage <- function(stage) {
        length(unique(as.character(
          membership$reaction_id[
            !is.na(membership$inclusion_stage) &
              membership$inclusion_stage == stage
          ]
        )))
      }

      answer$summary$n_core_reactions[[i]] <- length(unique(
        as.character(core_reactions$reaction_id[selected])
      ))
      answer$summary$n_core_genes[[i]] <- length(unique(
        as.character(core_reactions$gene[gene_selected])
      ))
      answer$summary$n_reactions[[i]] <- length(unique(
        as.character(membership$reaction_id)
      ))
      answer$summary$n_subsystem_added[[i]] <- count_stage(
        "same_core_subsystem"
      )
      answer$summary$n_database_added[[i]] <- count_stage(
        "shared_kegg_or_reactome_reaction"
      )
      answer$summary$n_rhea_added[[i]] <- count_stage(
        "shared_master_rhea_reaction"
      )
    }
  }

  answer
}
