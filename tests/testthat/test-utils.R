test_that("Rolling sum works", {
  n <- 100L
  test_data <- runif(n) * 10
  test_data[2] <- NA
  test_data[10] <- NA
  
  expand_grid(
    align = c("left", "center", "right"),
    k = 3:5,
    na_pad = c(TRUE, FALSE),
    na_rm = c(TRUE, FALSE)
  ) %>%
    pwalk(function(align, k, na_pad, na_rm) {
      if (!na_pad) {
        target <-
          zoo::rollsum(test_data,
                       k,
                       na.rm = na_rm,
                       align = align)
      } else {
        target <-
          zoo::rollsum(
            test_data,
            k,
            fill = ifelse(na_pad, NA, NULL),
            na.rm = na_rm,
            align = align
          )
      }
      result <- rollsum(test_data,
                        k,
                        na_pad = na_pad,
                        na.rm = na_rm,
                        align = align)
      # browser()
      expect_equal(result, target)
    })
})


test_that("Rolling sum works", {
  n <- 100L
  test_data <- runif(n) * 10
  test_data[2] <- NA
  test_data[10] <- NA
  
  expand_grid(
    align = c("left", "center", "right"),
    k = 3:5,
    na_pad = c(TRUE, FALSE),
    na_rm = c(TRUE, FALSE)
  ) %>%
    pwalk(function(align, k, na_pad, na_rm) {
      if (!na_pad) {
        target <-
          zoo::rollmean(test_data,
                        k,
                        na.rm = na_rm,
                        align = align)
      } else {
        target <-
          zoo::rollmean(
            test_data,
            k,
            fill = ifelse(na_pad, NA, NULL),
            na.rm = na_rm,
            align = align
          )
      }
      result <- rollmean(test_data,
                         k,
                         na_pad = na_pad,
                         na.rm = na_rm,
                         align = align)
      # browser()
      expect_equal(result, target)
    })
})



test_that("Obtaining Seqinfo works", {
  expect_equal(get_seqinfo(NULL), NULL)
  
  seqinfo <- get_seqinfo("GRCh37")
  expect_true(is(seqinfo, "Seqinfo"))
  
  seqinfo <- get_seqinfo("hs37-1kg")
  expect_true(is(seqinfo, "Seqinfo"))
  
  seqinfo <-
    get_seqinfo(
      "https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/1kg_v37.chrom.sizes"
    )
  expect_true(is(seqinfo, "Seqinfo"))
  expect_equal(GenomeInfoDb::genome(seqinfo) %>% unique, "1kg_v37")
  
  seqinfo <-
    get_seqinfo(
      "https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/1kg_v37.chrom.sizes",
      genome_name = "foo"
    )
  expect_true(is(seqinfo, "Seqinfo"))
  expect_equal(GenomeInfoDb::genome(seqinfo) %>% unique, "foo")
})


test_that("Making windows works", {
  # Genome: hs37-1kg
  windows <-
    make_windows(window_size = 500e3L,
                 genome = "hs37-1kg",
                 chrom = "22")
  expect_equal(width(windows) %>% head(-1) %>% unique(), 500e3L)
  expect_equal(
    GenomicRanges::seqinfo(windows) %>% GenomeInfoDb::genome() %>% unique(),
    "hs37-1kg"
  )
  
  # Genome: GRCh37
  windows <-
    make_windows(window_size = 500e3L,
                 genome = "GRCh37",
                 chrom = c("1", "2"))
  expect_equal(seqnames(windows) %>% as.character() %>% unique(), c("1", "2"))
  expect_equal(GenomicRanges::seqinfo(windows) %>% GenomeInfoDb::genome() %>% unique(),
               "GRCh37")
  
  # Genome: custom chrom.sizes
  tempbed <- tempfile(fileext = ".tsv")
  on.exit(unlink(tempbed), add = TRUE)
  data.table::data.table(chrom = c("chr1", "chr2"), 
                         size = c(640, 880)) %>% 
    data.table::fwrite(file = tempbed, sep = "\t", col.names = FALSE)
    
  windows <- make_windows(window_size = 100,
                          genome = tempbed)
  expect_equal(seqnames(windows) %>% as.character() %>% unique(), c("chr1", "chr2"))
  expect_equal(length(windows), 16)
})


test_that("Converting between data.frame and GenomicRanges works", {
  data0 <- read_bed("example2.bed", genome = "hs37-1kg")
  data1 <- as.bedtorch_table(data0)
  data2 <- as.GenomicRanges(data1)
  expect_equal(data0, data2)
  
  data0 <- read_bed("example2.bed")
  data1 <- as.bedtorch_table(data0)
  data2 <- as.GenomicRanges(data1)
  expect_equal(data0, data2)
  
  dt <- read_bed("example2.bed", use_gr = FALSE, genome = "hs37-1kg")
  dt_df <- as.data.frame(dt)
  expect_equal(dt, as.bedtorch_table(dt))
  expect_equal(dt, as.bedtorch_table(dt_df))

  gr <- read_bed("example2.bed", genome = "hs37-1kg")
  expect_equal(gr, as.GenomicRanges(dt))
  expect_equal(gr, as.GenomicRanges(dt_df))
})
