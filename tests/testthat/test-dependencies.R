test_that("requirements.tsv lists every package the engines and the code call", {
  req <- bioana:::.ba_tsv("requirements")$package
  eng <- unlist(strsplit(bioana:::.ba_tsv("engines")$requires, ",", fixed = TRUE))
  expect_true(all(trimws(eng[nzchar(eng)]) %in% req))
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir))
  code <- unlist(lapply(list.files(rdir, "\\.R$", full.names = TRUE), readLines))
  used <- unique(sub("::$", "", unlist(regmatches(code, gregexpr("\\b[A-Za-z][A-Za-z0-9.]+::", code)))))
  used <- c(used, unlist(regmatches(code, gregexpr('(?<=getExportedValue\\(")[A-Za-z0-9.]+', code, perl = TRUE))))
  base <- c("stats", "utils", "bioana", "tools", "methods", "grDevices", "graphics", "BiocManager", "remotes")
  expect_equal(setdiff(used, c(req, base)), character())
  d <- bioana_check_dependencies("minimal")
  expect_true(all(d$installed))
  expect_error(bioana_check_dependencies("nope"), "Unknown profile")
})

test_that("requirements.txt lists every package in requirements.tsv", {
  root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = FALSE)
  f <- file.path(root, "requirements.txt")
  skip_if_not(file.exists(f), "requirements.txt not in this build")
  txt <- readLines(f); pk <- sub("\\s.*$", "", txt[nzchar(txt) & !startsWith(txt, "#")])
  expect_setequal(pk, unique(bioana:::.ba_tsv("requirements")$package))
})
