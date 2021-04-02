test_that("Detecting gzipped files works", {
  expect_true(is_gzip("a.gz"))
  expect_true(is_gzip("a.bgz"))
  expect_false(is_gzip("a.bed"))
  #  Multipe
  expect_true(all(is_gzip(c("a.gz", "b.bed")) == c(TRUE, FALSE)))
})


test_that("Normalizing tabix rangees works", {
  expect_equal(normalize_tabix_range("2:1-2"), "2:1-2")
  expect_equal(normalize_tabix_range("2"), "2:1-2100000000")
  expect_equal(normalize_tabix_range(c("2:1-2", "3")), c("2:1-2", "3:1-2100000000"))
  expect_error(normalize_tabix_range(":1"))
})


test_that("Post-processing BED table works", {
  data <- data.table::data.table(V1 = 21, V2 = 1:10, V3 = 21:30, V4 = rnorm(10))
  data2 <- post_process_table(data)
  expect_equal(colnames(data2), c("chrom", "start", "end", "score"))
  
  data <- data.table::data.table(chrom = 21, s = 1:10, e = 21:30, V4 = "name")
  data2 <- post_process_table(data)
  expect_equal(colnames(data2), c("chrom", "start", "end", "V4"))
  
  data <- data.table::data.table(V1 = 21, V2 = 1:10, V3 = 21:30, V4 = "")
  data2 <- post_process_table(data)
  expect_equal(colnames(data2), c("chrom", "start", "end", "feature"))
  
  example_url <- "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
  data <- data.table::fread(example_url)
  data2 <- post_process_table(data)
  # Check the index
  expect_equal(attr(data2, which = "sorted"), c("chrom", "start", "end"))
  
  # Check whether chrom is factor
  expect_true(is.factor(data2$chrom))
  chrom_levels <- levels(data2$chrom)
  expect_equal(chrom_levels, str_sort(chrom_levels, numeric = TRUE))
  
  # Check start and end
  expect_true(data2[, all(list(start, end) %>% map_lgl(is.integer))])
})
