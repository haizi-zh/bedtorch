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

test_that("Reading data works", {
  data <- read_bed("example2.bed", use_gr = TRUE, genome = "GRCh37")
  expect_equal(length(data), 3000L)
  expect_equal(names(mcols(data)), c("score1", "score2"))
  expect_false(is.null(GenomeInfoDb::seqinfo(data)))
  
  expect_warning(data <- read_bed("example2.bed", use_gr = TRUE, genome = "GRCh37", range = "1:1101-1300"))
  expect_equal(length(data), 3L)
  expect_warning(data <- read_bed("example2.bed", use_gr = TRUE, genome = "GRCh37", range = "1:1100-1300"))
  expect_equal(length(data), 4L)
  expect_warning(data <- read_bed("example2.bed", use_gr = TRUE, genome = "GRCh37", range = "1:1100-1212"))
  expect_equal(length(data), 4L)
  expect_warning(data <- read_bed("example2.bed", use_gr = TRUE, genome = "GRCh37", range = "1:1100-1211"))
  expect_equal(length(data), 3L)
  
  data <- read_bed("example2.bed.gz", use_gr = TRUE, genome = "GRCh37", range = "1:1101-1300")
  expect_equal(names(mcols(data)), c("score1", "score2"))
  expect_equal(length(data), 3L)
})


test_that("Reading non-GenomicRanges data works", {
  data <- read_bed("example2.bed", use_gr = FALSE)
  expect_equal(colnames(data),
               c("chrom", "start", "end", "score1", "score2"))
  expect_equal(map_chr(data, class) %>% as.character(), c("factor", rep("integer", 4)))
})


test_that("Reading non-standard data works", {
  # Non-standard separator
  data <-
    read_bed(
      "https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt",
      sep = "auto",
      use_gr = FALSE
    )
  expect_equal(nrow(data), 26421)
  expect_equal(ncol(data), 5)
})


test_that("Reading remote data works", {
  url <- "https://github.com/haizi-zh/bedtorch/raw/4027bce5/inst/extdata/example2.bed.gz"
  index_url <- "https://github.com/haizi-zh/bedtorch/raw/4027bce5/inst/extdata/example2.bed.gz.tbi"
  
  data <- read_bed(url, use_gr = TRUE, genome = "GRCh37")
  expect_equal(length(data), 3000L)
  expect_equal(names(mcols(data)), c("score1", "score2"))
  expect_false(is.null(GenomeInfoDb::seqinfo(data)))
  
  data <- read_bed(url, tabix_index = index_url, use_gr = TRUE, genome = "GRCh37", range = "1:1101-1300")
  expect_equal(names(mcols(data)), c("score1", "score2"))
  expect_equal(length(data), 3L)
  
  url <- "https://raw.githubusercontent.com/haizi-zh/bedtorch/main/inst/extdata/example_merge.bed"
  data <- read_bed(url, use_gr = TRUE, genome = "GRCh37")
  expect_equal(length(data), 20L)
  
  # Non-standard URL
  url <- "https://git.io/JYATB"
  index_url <- "https://git.io/JYAkT"
  data <- read_bed(file_path =url)
  expect_equal(length(data), 5308L)
  
  data <- read_bed(file_path = url, tabix_index = index_url, range = "22:20000001-30000001")
  expect_equal(length(data), 21L)
})


test_that("Write data to files works", {
  dt <- read_bed("example2.bed")
  
  temp_bed1 <- tempfile(fileext = ".bed")
  on.exit(unlink(temp_bed1), add = TRUE)
  write_bed(dt, file_path = temp_bed1)
  expect_true(file.exists(temp_bed1))
  expect_equal(read_bed(temp_bed1), dt)
  
  temp_bed2 <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(temp_bed2), add = TRUE)
  write_bed(dt, file_path = temp_bed2, tabix_index = FALSE)
  expect_true(file.exists(temp_bed2))
  expect_false(file.exists(paste0(temp_bed2, ".tbi")))
  
  temp_bed3 <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(temp_bed3), add = TRUE)
  on.exit(unlink(paste0(temp_bed3, ".tbi")), add = TRUE)
  write_bed(dt, file_path = temp_bed3, tabix_index = TRUE)
  expect_true(file.exists(temp_bed3))
  expect_true(file.exists(paste0(temp_bed3, ".tbi")))
  
  temp_bed4 <- tempfile(fileext = ".bed.gz")
  on.exit(unlink(temp_bed4), add = TRUE)
  write_bed(dt, file_path = temp_bed4, tabix_index = FALSE, batch_size = 350)
  expect_equal(read_bed(temp_bed4), dt)
})