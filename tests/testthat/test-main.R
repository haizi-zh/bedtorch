library(here)

test_that("Loading a headless bedGraph file works", {
  data <- load_genbed(here("data-raw/example-01.bedGraph"), na = ".")
  expect_equal(nrow(data), 71L)
  expect_true(all(colnames(data) == c("chrom", "start", "end", "X4")))
  expect_true(all(summarize_all(data, class) == c("character", "integer", "integer", "numeric")))
  expect_equal(median(data$start), 33500000L)
  expect_equal(mean(data$X4, na.rm = TRUE), 0.05022998, tolerance = 1e-6)
})


test_that("Setting column names to a headless bedGraph file works", {
  names <- c("chr", "s", "e", "value")
  data <- load_genbed(here("data-raw/example-01.bedGraph"), na = ".", col_names = names)
  expect_equal(nrow(data), 71L)
  expect_true(all(colnames(data) == names))
  expect_true(all(summarize_all(data, class) == c("character", "integer", "integer", "numeric")))
  expect_equal(median(data$s), 33500000L)
  expect_equal(mean(data$value, na.rm = TRUE), 0.05022998, tolerance = 1e-6)
})


test_that("Loading bedGraph file with header works", {
  data <- load_genbed(here("data-raw/example-03.bedGraph"), na = ".")
  expect_equal(nrow(data), 71L)
  expect_true(all(colnames(data) == c("chr", "start", "end", "count")))
  expect_true(all(summarize_all(data, class) == c("character", "integer", "integer", "numeric")))
  expect_equal(median(data$start), 33500000L)
  expect_equal(mean(data$count, na.rm = TRUE), 0.05022998, tolerance = 1e-6)
})


test_that("Checking if a URL is remote works", {
  schema <- c("ftp", "sftp", "http", "https", "file")
  url <- c(paste0(schema, "://foo.com/bar"), "foo.bar")
  is_remote_url(url)
  expect_true(all(sapply(url, is_remote_url) == c(rep(TRUE, 4), rep(FALSE, 2))))
  expect_false(is_remote_url("foo.bar"))
})


test_that("Checking if a file is gzipped works", {
  expect_true(all(is_gzip_file(c("foo.bar", "foo.gz", "foo.bar.gz", "foo", "")) == c(FALSE, TRUE, TRUE, FALSE, FALSE)))
  expect_true(is.na(is_gzip_file(NA)))
})
