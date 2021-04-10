test_that("Merging overlapping intervals works", {
  dt <- read_bed(("example_merge.bed"))

  merged <- merge_bed(dt, max_dist = 2)
  expect_equal(merged, read_bed("example_merge_r2.bed"))

  merged <- merge_bed(dt)
  expect_true(all(merged == read_bed(("example_merge_r1.bed"))))

  merged <-
    merge_bed(dt, max_dist = 0, operation = list(score = list(on = "score", func = sum)))
  expect_equal(merged, read_bed("example_merge_r3.bed"))
  
  dt <- read_bed(("example2.bed"))
  
  merged <- merge_bed(dt)
  expect_equal(merged, read_bed(("example2_merge_r1.bed")))

  merged <- merge_bed(dt, max_dist = 5)
  expect_equal(merged, read_bed(("example2_merge_r2.bed")))
  
  merged <-
    merge_bed(
      dt,
      max_dist = 5,
      operation = list(score1 = list(on = "score1", func = function(x) sum(x, na.rm = TRUE)),
                       score2 = list(on = "score2", func = length))
    )
  target <- read_bed(("example2_merge_r3.bed"))
  names(mcols(target)) <- names(mcols(merged))
  expect_equal(merged, target)
})


test_that("Finding intervals overlapping with windows works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- intersect_bed(dt1, dt2, mode = "unique")
  target <- read_bed(("example2_intersect_r3.bed"))
  names(mcols(target)) <- names(mcols(result))
  expect_equal(result, target)
})


test_that("Finding intersections in default mode works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- sort(intersect_bed(dt1, dt2))
  target <- sort(read_bed("example2_intersect_r1.bed"))
  names(mcols(target)) <- names(mcols(result))
  expect_equal(result, target)
  
  # Set min_overlap
  result <- sort(intersect_bed(dt1, dt2, min_overlap = 10, min_overlap_type = "bp"))
  target <- read_bed("example2_intersect_r1.bed")
  target <- sort(target[width(target) >= 10])
  names(mcols(target)) <- names(mcols(result))
  expect_equal(result, target)
})


test_that("Finding intersections between two tables works", {
  # dt1 <- read_bed(("example_merge.bed")
  # dt2 <- read_bed(("example_intersect_y.bed")
  # 
  # result <- intersect_bed(dt1, dt2)
  # expect_true(all(result == read_bed(("example_intersect_r1.bed")))
  # 
  # result <- intersect_bed(dt1, dt2, mode = "unique")
  # expect_equal(result, read_bed(("example_intersect_r3.bed"))
  
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- intersect_bed(dt1, dt2)
  target <- read_bed(("example2_intersect_r1.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "unique")
  target <- read_bed(("example2_intersect_r3.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wa")
  target <- read_bed(("example2_intersect_r4_wa.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wb")
  target <- read_bed(("example2_intersect_r4_wb.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  target[, chrom.1 := factor(chrom.1, levels = levels(result$chrom.1))]
  setkey(result, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  setkey(target, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wo")
  target <- read_bed(("example2_intersect_r4_wo.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  target[, chrom.1 := factor(chrom.1, levels = levels(result$chrom.1))]
  setkey(result, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  setkey(target, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  expect_equal(result, target)
})

test_that("Finding intersections between two tables with min_overlap settings works", {
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- intersect_bed(dt1, dt2, min_overlap = 10, min_overlap_type = "bp")
  target <- read_bed(("example2_intersect_r1.bed"), use_gr = FALSE)[end - start >= 10]
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, min_overlap = 0.5, min_overlap_type = "frac1")
  target <- read_bed(("example2_intersect_r1_f0_5.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, min_overlap = 0.1, min_overlap_type = "frac2")
  target <- read_bed(("example2_intersect_r1_F0_1.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})

# test_that("Finding intersections between two tables with max_gap settings works", {
#   dt1 <- read_bed(("example_merge.bed")
#   dt2 <- read_bed(("example_intersect_y.bed")
#   
#   result <- intersect_bed(dt1, dt2, max_gap = 1)
#   expect_equal(nrow(result), nrow(read_bed(("example_intersect_r1.bed")) + 2)
#   
#   result_default <- intersect_bed(dt1, dt2)
#   result <- intersect_bed(dt1, dt2, min_overlap = 3)
#   expect_equal(result, result[end - start >= 3])
#   
#   dt1 <- read_bed(("example2.bed")
#   dt2 <- read_bed(("example2_window.bed")
#   
#   result <- intersect_bed(dt1, dt2, mode = "wo", max_gap = 10)
#   result <- result[]
#   target <- read_bed(("example2_intersect_r4_wa.bed")
#   setnames(target, new = colnames(result))
#   expect_equal(result, target)
# })

test_that("Excluding regions works (GenomicRanges)", {
  dt1 <- read_bed("example_merge.bed")
  dt2 <- read_bed("example_intersect_y.bed")
  
  result <- exclude_bed(dt1, dt2)
  expect_equal(result, read_bed("example_intersect_r2.bed"))
  
  dt1 <- read_bed("example2.bed")
  dt2 <- read_bed("example2_window.bed")
  
  result <- exclude_bed(dt1, dt2)
  target <- read_bed("example2_intersect_r2.bed")
  names(mcols(target)) <- names(mcols(result))
  expect_equal(result, target)
})

test_that("Excluding regions works", {
  dt1 <- read_bed(("example_merge.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example_intersect_y.bed"), use_gr = FALSE)
  
  result <- intersect_bed(dt1, dt2, mode = "exclude")
  expect_equal(result, read_bed(("example_intersect_r2.bed"), use_gr = FALSE))
  
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- intersect_bed(dt1, dt2, mode = "exclude")
  target <- read_bed(("example2_intersect_r2.bed"), use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})

test_that("Slopping intervals works", {
  dt <- read_bed("example_slop.bed", use_gr = FALSE)
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "both")
  target <- read_bed("example_slop_r1.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "right")
  target <- read_bed("example_slop_r2.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "left")
  target <- read_bed("example_slop_r3.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- slop_bed(dt, slop = 1.2, slop_type = "frac", mode = "both")
  target <- read_bed("example_slop_r4.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})


test_that("Mapping intervals works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- map_bed(dt1, dt2, operation = list(score1 = list(on = "score1", func = sum),
                                               score2 = list(on = "score2", func = length)))
  
  target <- read_bed(("example2_map_r1.bed"))
  target <- target[!is.na(target$score)]
  names(mcols(target)) <- names(mcols(result))
  
  expect_equal(result, target)
})


test_that("Subtracting BED tables works", {
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- subtract_bed(dt1, dt2)
  target <- read_bed("example2_subtract_r1.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})

test_that("Generating complement BED works", {
  dt1 <- read_bed(("example_merge.bed"), use_gr = FALSE)
  
  result <- complement_bed(dt1)
  target <- read_bed("example_complement_r1.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  
  result <- complement_bed(dt1)
  target <- read_bed("example_complement_r2.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})


test_that("Shuffling BED works", {
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- shuffle_bed(dt1, excluded_region = dt2)
  # Cross-chrom check
  expect_equal(dt1[, .N, by = chrom][[2]], result[, .N, by = chrom][[2]])
  expect_equal(dt1[, sum(end - start), by = chrom][[2]], result[, sum(end - start), by = chrom][[2]])
  # Excluded region check
  expect_equal(nrow(intersect(result, dt2)), 0)
  
  # Without excluded_region
  result <- shuffle_bed(dt1[1:100])
  # Cross-chrom check
  expect_equal(dt1[1:100, .N, by = chrom][[2]], result[, .N, by = chrom][[2]])
  expect_equal(dt1[1:100, sum(end - start), by = chrom][[2]], result[, sum(end - start), by = chrom][[2]])
  
  # Set RNG seed
  
  expect_equal(shuffle_bed(dt1[1:100], seed = 1),
               shuffle_bed(dt1[1:100], seed = 1))
})

test_that("Clustering BED intervals works", {
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  
  result <- cluster_bed(dt1)
  target <- read_bed("example2_cluster_r1.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- cluster_bed(dt1, max_dist = 10)
  target <- read_bed("example2_cluster_r2.bed", use_gr = FALSE)
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})


test_that("Calculating Jaccard distance works", {
  dt1 <- read_bed(("example2.bed"), use_gr = FALSE)
  dt2 <- read_bed(("example2_window.bed"), use_gr = FALSE)
  
  result <- jaccard_bed(dt1, dt2)
  expect_equal(result$intersection[1], 25241)
  expect_equal(result$union[1], 227080)
  expect_equal(result$n_intersections[1], 141)
  expect_true(abs(result$jaccard[1] - 0.111155) < 1e-5)
})


test_that("Making windows over genome works", {
  window_size <- 500e3L
  dt <- make_windows(window_size, "GRCh37", chrom = c("21", "22"))
  
  chrom_list <- as.character(unique(seqnames(dt)))
  expect_true(all(chrom_list == c("21", "22")))
  
  # Check the interval widths
  chrom_list %>% map_lgl(function(chrom) {
    dt <- dt[seqnames(dt) == chrom]
    all(head(width(dt), -1) == window_size)
  }) %>% all() %>% expect_true()
})
