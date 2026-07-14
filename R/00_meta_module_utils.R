.rc_mm_trim_unique <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}

.rc_mm_split_values <- function(x) {
  x <- as.character(x)
  out <- unlist(strsplit(x[!is.na(x) & nzchar(x)], "[;,|]", perl = TRUE), use.names = FALSE)
  .rc_mm_trim_unique(out)
}

.rc_mm_first_column <- function(x, candidates) {
  hit <- intersect(candidates, colnames(x))
  if (length(hit)) hit[[1L]] else NULL
}

.rc_mm_write_tsv_gz <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(file, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
  invisible(file)
}

.rc_mm_empty_edges <- function() {
  data.frame(sample_id = character(), gene_a = character(), gene_b = character(),
             edge_type = character(), shared_tf_count = integer(),
             projection_weight = numeric(), tf_jaccard = numeric(),
             direct_regulatory = logical(), module_id = character(),
             stringsAsFactors = FALSE)
}

.rc_mm_components <- function(nodes, edges) {
  nodes <- .rc_mm_trim_unique(nodes)
  if (!length(nodes)) return(data.frame(gene = character(), component = integer(), stringsAsFactors = FALSE))
  parent <- stats::setNames(nodes, nodes)
  find_root <- function(x) {
    y <- x
    while (!identical(parent[[y]], y)) y <- parent[[y]]
    y
  }
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      a <- as.character(edges$gene_a[[i]])
      b <- as.character(edges$gene_b[[i]])
      if (!a %in% nodes || !b %in% nodes) next
      ra <- find_root(a)
      rb <- find_root(b)
      if (!identical(ra, rb)) parent[[rb]] <- ra
    }
  }
  roots <- vapply(nodes, find_root, character(1))
  root_levels <- unique(roots)
  data.frame(gene = nodes, component = match(roots, root_levels), stringsAsFactors = FALSE)
}


.rc_validate_pando_install <- function(pando_version = NULL,
                                       pando_remote_username = "1667857557",
                                       pando_remote_repo = "Pando_regcompass",
                                       require_pando_remote = TRUE,
                                       description = NULL,
                                       installed_version = NULL) {
  if (is.null(installed_version)) installed_version <- as.character(utils::packageVersion("Pando"))
  installed_version <- as.character(installed_version)
  if (is.null(description)) description <- utils::packageDescription("Pando")
  required_user <- .rc_mm_trim_unique(pando_remote_username)
  required_repo <- .rc_mm_trim_unique(pando_remote_repo)
  remote_username <- as.character(description$RemoteUsername %||% description$GithubUsername %||% NA_character_)
  remote_repo <- as.character(description$RemoteRepo %||% description$GithubRepo %||% NA_character_)
  remote_ref <- as.character(description$RemoteRef %||% description$GithubRef %||% NA_character_)
  remote_sha <- as.character(description$RemoteSha %||% description$GithubSHA1 %||% NA_character_)
  if (isTRUE(require_pando_remote)) {
    if (!length(required_user) || is.na(remote_username) || !remote_username %in% required_user) {
      stop("Pando remote username mismatch: installed ", remote_username,
           ", required one of ", paste(required_user, collapse = ", "), ".", call. = FALSE)
    }
    if (!length(required_repo) || is.na(remote_repo) || !remote_repo %in% required_repo) {
      stop("Pando remote repository mismatch: installed ", remote_repo,
           ", required one of ", paste(required_repo, collapse = ", "), ".", call. = FALSE)
    }
  }
  list(version = installed_version,
       remote_username = remote_username,
       remote_repo = remote_repo,
       remote_ref = remote_ref,
       remote_sha = remote_sha)
}

