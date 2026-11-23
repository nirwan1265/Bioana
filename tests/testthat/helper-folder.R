# A copy of the 2-treatment x 6-block RCBD fixture (a Biodes design folder) to fill with data.
rcbd_folder <- function() {
  d <- file.path(tempfile("bioana_"), "design")
  dir.create(d, recursive = TRUE)
  file.copy(list.files(testthat::test_path("fixtures", "rcbd"), full.names = TRUE), d)
  d
}
read_meta <- function(d) utils::read.delim(file.path(d, "sample_metadata_template.tsv"), colClasses = "character", check.names = FALSE)
write_meta <- function(m, d) utils::write.table(m, file.path(d, "sample_metadata_template.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write_mat <- function(m, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(data.frame(feature_id = rownames(m), m, check.names = FALSE), path, sep = "\t", row.names = FALSE, quote = FALSE)
}
