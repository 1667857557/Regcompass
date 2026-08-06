args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args

find_root <- function() {
  workspace <- Sys.getenv("GITHUB_WORKSPACE", unset = "")
  candidates <- unique(c(
    if (nzchar(workspace)) workspace else character(),
    ".", "..", file.path("..", "..")
  ))
  candidates <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (!length(candidates)) stop("Could not locate repository root.", call. = FALSE)
  normalizePath(candidates[[1L]], mustWork = TRUE)
}

root <- find_root()
old_wd <- setwd(root)
on.exit(setwd(old_wd), add = TRUE)

description <- read.dcf("DESCRIPTION")
collate_text <- unname(description[1L, "Collate"])
collate <- regmatches(
  collate_text,
  gregexpr("'[^']+[.]R'", collate_text, perl = TRUE)
)[[1L]]
collate <- gsub("^'|'$", "", collate)
collate_paths <- file.path("R", collate)
missing_collate <- collate_paths[!file.exists(collate_paths)]
all_r_files <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
uncollated <- setdiff(all_r_files, collate_paths)

namespace <- readLines("NAMESPACE", warn = FALSE)
exports <- sub(
  "^export[(]([^)]*)[)].*$", "\\1",
  grep("^export[(]", namespace, value = TRUE)
)
s3_methods <- sub(
  "^S3method[(]([^,]+),([^)]*)[)].*$", "\\1.\\2",
  grep("^S3method[(]", namespace, value = TRUE)
)

extract_definition <- function(expr, file) {
  if (!is.call(expr) || length(expr) < 3L) return(NULL)
  op <- as.character(expr[[1L]])
  if (!op %in% c("<-", "=")) return(NULL)
  lhs <- expr[[2L]]
  rhs <- expr[[3L]]
  if (!is.symbol(lhs) || !is.call(rhs) || !identical(rhs[[1L]], as.name("function"))) {
    return(NULL)
  }
  ref <- attr(expr, "srcref")
  line <- if (!is.null(ref)) as.integer(ref[[1L]]) else NA_integer_
  body_expr <- rhs[[3L]]
  if (is.call(body_expr) && identical(body_expr[[1L]], as.name("{")) &&
      length(body_expr) == 2L) {
    body_expr <- body_expr[[2L]]
  }
  forward_target <- NA_character_
  if (is.call(body_expr)) {
    head <- body_expr[[1L]]
    if (is.symbol(head)) {
      candidate <- as.character(head)
      if (!candidate %in% c(
        "if", "for", "while", "repeat", "return", "stop", "warning",
        "tryCatch", "withCallingHandlers", "switch", "{", "(", "<-", "="
      )) forward_target <- candidate
    } else if (is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
      forward_target <- paste0(
        as.character(head[[2L]]), as.character(head[[1L]]),
        as.character(head[[3L]])
      )
    }
  }
  data.frame(
    name = as.character(lhs),
    file = file,
    line = line,
    forward_target = forward_target,
    stringsAsFactors = FALSE
  )
}

definitions <- list()
for (file in all_r_files) {
  parsed <- tryCatch(
    parse(file = file, keep.source = TRUE),
    error = function(e) stop("Failed to parse ", file, ": ", conditionMessage(e),
                             call. = FALSE)
  )
  for (expr in parsed) {
    one <- extract_definition(expr, file)
    if (!is.null(one)) definitions[[length(definitions) + 1L]] <- one
  }
}
definitions <- if (length(definitions)) do.call(rbind, definitions) else {
  data.frame(
    name = character(), file = character(), line = integer(),
    forward_target = character(), stringsAsFactors = FALSE
  )
}
rownames(definitions) <- NULL

count_fixed <- function(name, files) {
  files <- files[file.exists(files)]
  if (!length(files)) return(0L)
  total <- 0L
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    matches <- gregexpr(name, lines, fixed = TRUE)
    total <- total + sum(vapply(matches, function(x) {
      if (length(x) == 1L && x[[1L]] == -1L) 0L else length(x)
    }, integer(1)))
  }
  total
}

test_files <- sort(unique(c(
  list.files("tests", pattern = "[.](R|Rmd)$", recursive = TRUE,
             full.names = TRUE),
  list.files("vignettes", pattern = "[.](R|Rmd)$", recursive = TRUE,
             full.names = TRUE)
)))
doc_files <- sort(unique(c(
  "README.md", "NEWS.md", "DESCRIPTION", "NAMESPACE",
  list.files("docs", pattern = "[.](md|Rmd)$", recursive = TRUE,
             full.names = TRUE),
  list.files("man", pattern = "[.]Rd$", recursive = TRUE,
             full.names = TRUE),
  list.files(".github", pattern = "[.](yml|yaml|md)$", recursive = TRUE,
             full.names = TRUE)
)))

if (nrow(definitions)) {
  definitions$definition_count <- ave(
    definitions$name, definitions$name, FUN = length
  )
  unique_defs <- definitions[!duplicated(definitions$name), , drop = FALSE]
  unique_defs$source_occurrences <- vapply(
    unique_defs$name, count_fixed, integer(1), files = all_r_files
  )
  unique_defs$test_occurrences <- vapply(
    unique_defs$name, count_fixed, integer(1), files = test_files
  )
  unique_defs$doc_occurrences <- vapply(
    unique_defs$name, count_fixed, integer(1), files = doc_files
  )
  unique_defs$exported <- unique_defs$name %in% exports
  unique_defs$s3_registered <- unique_defs$name %in% s3_methods
  unique_defs$compatibility_named <- grepl(
    "compat|legacy|deprecated|obsolete|old", unique_defs$name,
    ignore.case = TRUE
  ) | grepl(
    "compat|legacy|deprecated|obsolete|old", unique_defs$file,
    ignore.case = TRUE
  )
  special <- c(
    ".onLoad", ".onUnload", ".onAttach", ".Last.lib", ".First.lib"
  )
  unique_defs$private_zero_reference <-
    !unique_defs$exported & !unique_defs$s3_registered &
    !unique_defs$name %in% special &
    unique_defs$source_occurrences <= unique_defs$definition_count &
    unique_defs$test_occurrences == 0L &
    unique_defs$doc_occurrences == 0L
  unique_defs$exported_without_usage <-
    unique_defs$exported &
    unique_defs$source_occurrences <= unique_defs$definition_count &
    unique_defs$test_occurrences == 0L &
    unique_defs$doc_occurrences <= 2L
} else {
  unique_defs <- definitions
}

duplicates <- definitions[duplicated(definitions$name) |
  duplicated(definitions$name, fromLast = TRUE), , drop = FALSE]
undefined_exports <- setdiff(exports, definitions$name)
undocumented_exports <- exports[!vapply(exports, function(name) {
  any(vapply(list.files("man", pattern = "[.]Rd$", full.names = TRUE),
             function(file) {
               any(grepl(paste0("\\\\alias\\{", name, "\\}"),
                         readLines(file, warn = FALSE)))
             }, logical(1)))
}, logical(1))]

marker_files <- c(all_r_files, test_files, doc_files)
marker_patterns <- c(
  "[.]Deprecated[(]", "[.]Defunct[(]", "deprecated", "obsolete",
  "legacy API", "compatibility wrapper"
)
marker_hits <- list()
for (pattern in marker_patterns) {
  for (file in marker_files[file.exists(marker_files)]) {
    lines <- readLines(file, warn = FALSE)
    hit <- grep(pattern, lines, ignore.case = TRUE, perl = TRUE)
    if (length(hit)) {
      marker_hits[[length(marker_hits) + 1L]] <- data.frame(
        pattern = pattern, file = file, line = hit,
        text = trimws(lines[hit]), stringsAsFactors = FALSE
      )
    }
  }
}
marker_hits <- if (length(marker_hits)) do.call(rbind, marker_hits) else {
  data.frame(
    pattern = character(), file = character(), line = integer(),
    text = character(), stringsAsFactors = FALSE
  )
}

cat("API surface audit\n")
cat("=================\n")
cat("R files:", length(all_r_files), "\n")
cat("Top-level functions:", nrow(unique_defs), "\n")
cat("Exports:", length(exports), "\n\n")

print_section <- function(title, value) {
  cat("\n", title, "\n", strrep("-", nchar(title)), "\n", sep = "")
  if (is.data.frame(value)) {
    if (!nrow(value)) cat("none\n") else print(value, row.names = FALSE)
  } else {
    if (!length(value)) cat("none\n") else cat(paste(value, collapse = "\n"), "\n")
  }
}

print_section("Missing Collate files", missing_collate)
print_section("Uncollated R files", uncollated)
print_section("Duplicate top-level function definitions", duplicates)
print_section("Undefined exports", undefined_exports)
print_section("Undocumented exports", undocumented_exports)
print_section(
  "Private functions with zero repository references",
  unique_defs[unique_defs$private_zero_reference %in% TRUE,
              c("name", "file", "line", "forward_target",
                "source_occurrences", "test_occurrences", "doc_occurrences"),
              drop = FALSE]
)
print_section(
  "Compatibility-named functions/files",
  unique_defs[unique_defs$compatibility_named %in% TRUE,
              c("name", "file", "line", "forward_target",
                "source_occurrences", "test_occurrences", "doc_occurrences",
                "exported"), drop = FALSE]
)
print_section(
  "Single-call forwarding functions",
  unique_defs[!is.na(unique_defs$forward_target),
              c("name", "file", "line", "forward_target",
                "source_occurrences", "test_occurrences", "doc_occurrences",
                "exported"), drop = FALSE]
)
print_section(
  "Exported functions without non-man usage",
  unique_defs[unique_defs$exported_without_usage %in% TRUE,
              c("name", "file", "line", "source_occurrences",
                "test_occurrences", "doc_occurrences"), drop = FALSE]
)
print_section("Deprecation/obsolete markers", marker_hits)

hard_fail <- length(missing_collate) || length(uncollated) ||
  nrow(duplicates) || length(undefined_exports) || length(undocumented_exports) ||
  any(unique_defs$private_zero_reference %in% TRUE)

if (strict && hard_fail) {
  stop("API surface audit found hard failures; inspect the sections above.",
       call. = FALSE)
}

cat("\nAPI surface audit completed.\n")
