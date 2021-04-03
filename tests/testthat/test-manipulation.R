test_that("Merging overlapping intervals works", {
  dt <- read_bed(("example_merge.bed"))

  merged <- merge_bed(dt, max_dist = 2)
  expect_true(all(merged == read_bed(("example_merge_r2.bed"))))

  merged <- merge_bed(dt)
  expect_true(all(merged == read_bed(("example_merge_r1.bed"))))

  merged <-
    merge_bed(dt, max_dist = 0, operation = list(score = function(v) sum(v$score)))
  expect_equal(merged, read_bed(("example_merge_r3.bed")))
  
  dt <- read_bed(("example2.bed"))
  
  merged <- merge_bed(dt)
  expect_equal(merged, read_bed(("example2_merge_r1.bed")))

  merged <- merge_bed(dt, max_dist = 5)
  expect_equal(merged, read_bed(("example2_merge_r2.bed")))
  
  merged <-
    merge_bed(
      dt,
      max_dist = 5,
      operation = list(
        score1 = function(v)
          as.numeric(sum(v$score1, na.rm = TRUE)),
        score2 = function(v)
          as.numeric(length(v$score2))
      )
    )
  target <- read_bed(("example2_merge_r3.bed"))
  setnames(target, new = colnames(merged))
  expect_equal(merged, target)
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
  
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- intersect_bed(dt1, dt2)
  target <- read_bed(("example2_intersect_r1.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "unique")
  target <- read_bed(("example2_intersect_r3.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wa")
  target <- read_bed(("example2_intersect_r4_wa.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wb")
  target <- read_bed(("example2_intersect_r4_wb.bed"))
  setnames(target, new = colnames(result))
  target[, chrom.1 := factor(chrom.1, levels = levels(result$chrom.1))]
  setkey(result, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  setkey(target, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, mode = "wo")
  target <- read_bed(("example2_intersect_r4_wo.bed"))
  setnames(target, new = colnames(result))
  target[, chrom.1 := factor(chrom.1, levels = levels(result$chrom.1))]
  setkey(result, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  setkey(target, "chrom", "start", "end", "chrom.1", "start.1", "end.1")
  expect_equal(result, target)
})

test_that("Finding intersections between two tables with min_overlap settings works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- intersect_bed(dt1, dt2, min_overlap = 10, min_overlap_type = "bp")
  target <- read_bed(("example2_intersect_r1.bed"))[end - start >= 10]
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, min_overlap = 0.5, min_overlap_type = "frac1")
  target <- read_bed(("example2_intersect_r1_f0_5.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- intersect_bed(dt1, dt2, min_overlap = 0.1, min_overlap_type = "frac2")
  target <- read_bed(("example2_intersect_r1_F0_1.bed"))
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

test_that("Excluding regions works", {
  dt1 <- read_bed(("example_merge.bed"))
  dt2 <- read_bed(("example_intersect_y.bed"))
  
  result <- intersect_bed(dt1, dt2, mode = "exclude")
  expect_equal(result, read_bed(("example_intersect_r2.bed")))
  
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- intersect_bed(dt1, dt2, mode = "exclude")
  target <- read_bed(("example2_intersect_r2.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})

test_that("Slopping intervals works", {
  dt <- read_bed("example_slop.bed")
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "both")
  target <- read_bed("example_slop_r1.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "right")
  target <- read_bed("example_slop_r2.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- slop_bed(dt, slop = 10L, slop_type = "bp", mode = "left")
  target <- read_bed("example_slop_r3.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  result <- slop_bed(dt, slop = 1.2, slop_type = "frac", mode = "both")
  target <- read_bed("example_slop_r4.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})


test_that("Mapping intervals works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <-
    map_bed(dt1,
            dt2,
            operation = list(
              score1 = function(x)
                sum(x$score1),
              score2 = function(x)
                length(x$score2)
            ))
  result[is.na(score2), score2 := 0]
  target <- read_bed(("example2_map_r1.bed"))
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})


test_that("Subtracting BED tables works", {
  dt1 <- read_bed(("example2.bed"))
  dt2 <- read_bed(("example2_window.bed"))
  
  result <- subtract_bed(dt1, dt2)
  target <- read_bed("example2_subtract_r1.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})

test_that("Generating complement BED works", {
  dt1 <- read_bed(("example_merge.bed"))
  
  result <- complement_bed(dt1)
  target <- read_bed("example_complement_r1.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
  
  dt1 <- read_bed(("example2.bed"))
  
  result <- complement_bed(dt1)
  target <- read_bed("example_complement_r2.bed")
  setnames(target, new = colnames(result))
  expect_equal(result, target)
})
