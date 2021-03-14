library(here)

test_that("Loading a headless bedGraph file works", {
  for (file_path in c(
    here("data-raw/example-01.bedGraph"),
    "https://github.com/haizi-zh/bioessentials/raw/main/data-raw/example-01.bedGraph"
  )) {
    data <- load_genbed(file_path = file_path, na = ".")
    expect_equal(nrow(data), 71L)
    expect_true(all(colnames(data) == c("chrom", "start", "end", "X4")))
    expect_true(all(
      summarize_all(data, class) == c("character", "integer", "integer", "numeric")
    ))
    expect_equal(median(data$start), 33500000L)
    expect_equal(mean(data$X4, na.rm = TRUE), 0.05022998, tolerance = 1e-6)
  }
})


test_that("Loading a gzippped bedGraph file works", {
  for (file_path in c(
    here("data-raw/example-02.bedGraph.gz"),
    "https://github.com/haizi-zh/bioessentials/raw/main/data-raw/example-02.bedGraph.gz"
  )) {
    data <- load_genbed(file_path = file_path)
    expect_equal(nrow(data), 5308L)
    expect_true(all(
      colnames(data) == c("chrom", "start", "end", "score", "comp_method")
    ))
    expect_true(all(
      summarize_all(data, class) == c("character", "integer", "integer", "numeric", "character")
    ))
    expect_equal(mean(data$start), 80616334L)
    expect_equal(mean(data$score, na.rm = TRUE),-0.000330179264, tolerance = 1e-8)
  }
})


test_that("Selectively loading a gzippped bedGraph file works", {
  for (file_path in c(
    here("data-raw/example-02.bedGraph.gz"),
    "https://github.com/haizi-zh/bioessentials/raw/main/data-raw/example-02.bedGraph.gz"
  )) {
    data <-
      load_genbed(
        file_path = file_path,
        region = c("1:1000000-1000001", "3"),
        is_bedgraph = TRUE
      )
    expect_equal(nrow(data), 388)
    expect_true(all(
      colnames(data) == c("chrom", "start", "end", "score", "comp_method")
    ))
    expect_true(all(
      summarize_all(data, class) == c("character", "integer", "integer", "numeric", "character")
    ))
    expect_equal(mean(data$start), 98389175L)
    expect_equal(mean(data$score, na.rm = TRUE), 0.00069246274, tolerance = 1e-8)
  }
})


test_that("Selectively loading bedGraph file with header and NA values works",
          {
            data <-
              load_genbed(
                "https://github.com/haizi-zh/bioessentials/raw/main/data-raw/example-03.bedGraph.gz",
                region = "22:30000001-50000000",
                na = ".",
                is_bedgraph = TRUE
              )
            expect_equal(nrow(data), 40L)
            expect_true(all(colnames(data) == c("chr", "start", "end", "count")))
            expect_true(all(
              summarize_all(data, class) == c("character", "integer", "integer", "numeric")
            ))
            expect_equal(mean(data$start), 39750000L)
            expect_equal(mean(data$count, na.rm = TRUE), 0.04491535, tolerance = 1e-6)
          })


test_that("Setting column names to a headless bedGraph file works", {
  names <- c("chr", "s", "e", "value")
  data <-
    load_genbed(
      here("data-raw/example-01.bedGraph"),
      na = ".",
      col_names = names,
      is_bedgraph = TRUE
    )
  expect_equal(nrow(data), 71L)
  expect_true(all(colnames(data) == names))
  expect_true(all(
    summarize_all(data, class) == c("character", "integer", "integer", "numeric")
  ))
  expect_equal(median(data$s), 33500000L)
  expect_equal(mean(data$value, na.rm = TRUE), 0.05022998, tolerance = 1e-6)
})


test_that("Loading bedGraph file with header works", {
  data <-
    load_genbed(here("data-raw/example-03.bedGraph"),
                na = ".",
                is_bedgraph = TRUE)
  expect_equal(nrow(data), 71L)
  expect_true(all(colnames(data) == c("chr", "start", "end", "count")))
  expect_true(all(
    summarize_all(data, class) == c("character", "integer", "integer", "numeric")
  ))
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
  expect_true(all(is_gzip_file(
    c("foo.bar", "foo.gz", "foo.bar.gz", "foo", "")
  ) == c(FALSE, TRUE, TRUE, FALSE, FALSE)))
  expect_true(is.na(is_gzip_file(NA)))
})
