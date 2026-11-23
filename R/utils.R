`%||%` <- function(x, y) if (is.null(x)) y else x

# A tab-separated rules table from inst/rules (every column character).
.ba_tsv <- function(name) {
  p <- system.file("rules", paste0(name, ".tsv"), package = "bioana")
  if (!nzchar(p)) p <- file.path("inst", "rules", paste0(name, ".tsv"))
  utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(), quote = "")
}

.ba_write_tsv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  path
}

# Read a feature x sample table (first column = feature id), gz or not.
.ba_read_matrix <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "double"
  rownames(m) <- as.character(x[[1]])
  m
}
