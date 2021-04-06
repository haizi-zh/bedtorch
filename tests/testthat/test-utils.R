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