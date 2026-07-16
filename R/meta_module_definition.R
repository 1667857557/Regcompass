# GRN-to-reaction biological meta-module definition.

#' Map GRN metabolic gene modules to core Human-GEM reactions
rc_map_meta_module_core_reactions <- function(gene_nodes, gpr_table) {
  if (!is.data.frame(gene_nodes) ||
      !all(c("sample_id", "gene", "module_id") %in% colnames(gene_nodes))) {
    stop("`gene_nodes` must contain sample_id, gene and module_id.",
         call. = FALSE)
  }
  if (!is.data.frame(gpr_table) ||
      !all(c("reaction_id", "gene") %in% colnames(gpr_table))) {
    stop("`gpr_table` must contain reaction_id and gene.", call. = FALSE)
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
  if (!"and_group_id" %in% colnames(gpr)) {
    gpr$and_group_id <- ave(
      seq_len(nrow(gpr)),
      gpr$reaction_id,
      FUN = seq_along
    )
  }
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

#' Expand core reactions into GRN-defined biological reaction meta-modules
#'
#' Expansion is ordered: core subsystems, shared KEGG/Reactome identifiers, then
#' shared master-Rhea identifiers. The output is biological membership, not a
#' flux-feasibility support set.
.rc_expand_meta_module_reactions_core <- function(gem, core_reactions,
                                             subsystem_table = NULL,
                                             expansion_mode = c(
                                               "ordered_once", "fixed_point"
                                             ),
                                             max_iterations = 10L) {
  expansion_mode <- match.arg(expansion_mode)
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
  maps <- rc_reaction_crossref_maps(
    gem,
    subsystem_table = subsystem_table
  )
  maps$subsystem <- maps$subsystem[
  !is.na(maps$subsystem$reaction_id) &
    nzchar(trimws(as.character(maps$subsystem$reaction_id))) &
    !is.na(maps$subsystem$subsystem_id) &
    nzchar(trimws(as.character(maps$subsystem$subsystem_id))),
  , drop = FALSE
]
maps$kegg <- maps$kegg[
  !is.na(maps$kegg$reaction_id) &
    nzchar(trimws(as.character(maps$kegg$reaction_id))) &
    !is.na(maps$kegg$kegg_id) &
    nzchar(trimws(as.character(maps$kegg$kegg_id))),
  , drop = FALSE
]
maps$reactome <- maps$reactome[
  !is.na(maps$reactome$reaction_id) &
    nzchar(trimws(as.character(maps$reactome$reaction_id))) &
    !is.na(maps$reactome$reactome_id) &
    nzchar(trimws(as.character(maps$reactome$reactome_id))),
  , drop = FALSE
]
maps$rhea_master <- maps$rhea_master[
  !is.na(maps$rhea_master$reaction_id) &
    nzchar(trimws(as.character(maps$rhea_master$reaction_id))) &
    !is.na(maps$rhea_master$rhea_master_id) &
    nzchar(trimws(as.character(maps$rhea_master$rhea_master_id))),
  , drop = FALSE
]
  if (!nrow(maps$subsystem)) {
    stop(
      "No usable reaction-to-subsystem annotations were found.",
      call. = FALSE
    )
  }
  valid_reactions <- colnames(rc_validate_gem(gem)$S)
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
    stop(
      "No GRN genes mapped to valid GEM reactions.",
      call. = FALSE
    )
  }

  groups <- split(
    seq_len(nrow(core_reactions)),
    paste(
      core_reactions$sample_id,
      core_reactions$module_id,
      sep = "\001"
    )
  )
  membership_rows <- list()
  summary_rows <- list()

  for (key in names(groups)) {
    index <- groups[[key]]
    core <- unique(
      as.character(core_reactions$reaction_id[index])
    )
    sample_id <- as.character(
      core_reactions$sample_id[index[[1L]]]
    )
    module_id <- as.character(
      core_reactions$module_id[index[[1L]]]
    )
    included <- core
    reasons <- stats::setNames(
      rep("core_grn_gene", length(core)),
      core
    )
    source_ids <- stats::setNames(
      rep(NA_character_, length(core)),
      core
    )

    add_reactions <- function(reactions, reason, source_map = NULL) {
      reactions <- intersect(
        .rc_mm_trim_unique(reactions),
        valid_reactions
      )
      new <- setdiff(reactions, included)
      if (length(new)) {
        included <<- c(included, new)
        reasons[new] <<- reason
        if (!is.null(source_map)) {
          source_ids[new] <<- source_map[new]
        }
      }
      invisible(new)
    }

    expand_once <- function() {
      before <- length(included)
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
      source_for_subsystem_members <- function(
        subsystems,
        labels_by_subsystem
      ) {
        reactions <- reactions_in_subsystems(subsystems)
        stats::setNames(vapply(reactions, function(reaction) {
          subsystem_ids <- unique(maps$subsystem$subsystem_id[
            maps$subsystem$reaction_id == reaction &
              maps$subsystem$subsystem_id %in% subsystems
          ])
          labels <- unique(unlist(
            labels_by_subsystem[subsystem_ids],
            use.names = FALSE
          ))
          paste(
            labels[!is.na(labels) & nzchar(labels)],
            collapse = ";"
          )
        }, character(1)), reactions)
      }

      core_subsystems <- subsystems_for_reactions(core)
      subsystem_reactions <- reactions_in_subsystems(
        core_subsystems
      )
      subsystem_labels <- stats::setNames(
        lapply(
          core_subsystems,
          function(value) paste0("SUBSYSTEM:", value)
        ),
        core_subsystems
      )
      add_reactions(
        subsystem_reactions,
        "same_core_subsystem",
        source_for_subsystem_members(
          core_subsystems,
          subsystem_labels
        )
      )

      current <- included
      kegg_ids <- unique(maps$kegg$kegg_id[
        maps$kegg$reaction_id %in% current
      ])
      reactome_ids <- unique(maps$reactome$reactome_id[
        maps$reactome$reaction_id %in% current
      ])
      all_subsystems <- unique(maps$subsystem$subsystem_id)
      database_labels <- lapply(all_subsystems, function(subsystem) {
        reactions <- reactions_in_subsystems(subsystem)
        kegg <- intersect(
          unique(maps$kegg$kegg_id[
            maps$kegg$reaction_id %in% reactions
          ]),
          kegg_ids
        )
        reactome <- intersect(
          unique(maps$reactome$reactome_id[
            maps$reactome$reaction_id %in% reactions
          ]),
          reactome_ids
        )
        unique(c(
          paste0("SUBSYSTEM:", subsystem),
          paste0("KEGG:", kegg),
          paste0("REACTOME:", reactome)
        ))
      })
      names(database_labels) <- all_subsystems
      database_subsystems <- all_subsystems[vapply(
        database_labels,
        function(value) {
          any(grepl("^(KEGG|REACTOME):", value))
        },
        logical(1)
      )]
      add_reactions(
        reactions_in_subsystems(database_subsystems),
        "shared_kegg_or_reactome_subsystem",
        source_for_subsystem_members(
          database_subsystems,
          database_labels
        )
      )

      current <- included
      master_ids <- unique(maps$rhea_master$rhea_master_id[
        maps$rhea_master$reaction_id %in% current
      ])
      anchor_reactions <- unique(maps$rhea_master$reaction_id[
        maps$rhea_master$rhea_master_id %in% master_ids
      ])
      rhea_subsystems <- subsystems_for_reactions(
        anchor_reactions
      )
      rhea_labels <- lapply(rhea_subsystems, function(subsystem) {
        reactions <- reactions_in_subsystems(subsystem)
        identifiers <- intersect(
          unique(maps$rhea_master$rhea_master_id[
            maps$rhea_master$reaction_id %in% reactions
          ]),
          master_ids
        )
        unique(c(
          paste0("SUBSYSTEM:", subsystem),
          paste0("RHEA_MASTER:", identifiers)
        ))
      })
      names(rhea_labels) <- rhea_subsystems
      add_reactions(
        reactions_in_subsystems(rhea_subsystems),
        "shared_master_rhea_subsystem",
        source_for_subsystem_members(
          rhea_subsystems,
          rhea_labels
        )
      )
      length(included) > before
    }

    changed <- expand_once()
    iteration <- 1L
    if (identical(expansion_mode, "fixed_point")) {
      while (isTRUE(changed) &&
             iteration < as.integer(max_iterations)) {
        iteration <- iteration + 1L
        changed <- expand_once()
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
          "shared_kegg_or_reactome_subsystem"
      ),
      n_rhea_added = sum(
        membership$inclusion_stage ==
          "shared_master_rhea_subsystem"
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

rc_expand_meta_module_reactions <- function(gem, core_reactions,
                                             subsystem_table = NULL,
                                             expansion_mode = c(
                                               "ordered_once", "fixed_point"
                                             ),
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
  partial_anchor <- membership_keys %in%
    setdiff(candidate_keys, hard_keys)
  if (any(partial_anchor)) {
    answer$reaction_membership <- answer$reaction_membership[
      !partial_anchor,
      , drop = FALSE
    ]
    membership_keys <- key(answer$reaction_membership)
  }

  if (is.data.frame(answer$summary) && nrow(answer$summary)) {
    for (i in seq_len(nrow(answer$summary))) {
      selected <- as.character(core_reactions$sample_id) ==
        as.character(answer$summary$sample_id[[i]]) &
        as.character(core_reactions$module_id) ==
        as.character(answer$summary$module_id[[i]]) &
        core_reactions$is_core %in% TRUE
      answer$summary$n_core_reactions[[i]] <- length(unique(
        as.character(core_reactions$reaction_id[selected])
      ))
      gene_selected <- selected
      if ("group_complete" %in% colnames(core_reactions)) {
        gene_selected <- gene_selected & core_reactions$group_complete %in% TRUE
      }
      answer$summary$n_core_genes[[i]] <- length(unique(
        as.character(core_reactions$gene[gene_selected])
      ))
    }
  }
  answer
}
