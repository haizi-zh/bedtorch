# Find a new column name which does not conflict with existing names
# 
# # Returns "chrom.1"
# .available_colname(dt, "chrom")
.available_colname <- function(dt, name) {
  fields <- colnames(dt)
  
  suffix <- 0
  candidate <- name
  while (TRUE) {
    if (!candidate %in% fields)
      return(candidate)
    
    suffix <- suffix + 1
    candidate <- paste0(name, ".", suffix)
  }
}


#' Calculate rolling sums
#' 
#' Compared with [zoo::rollsum()], this function is implemented in C and is one
#' magnitude faster.
#' @param x A numeric vector.
#' @param k integer width of the rolling window.
#' @param na_pad Logical value about whether to pad left/right ends with NAs.
#' @param na.rm Similar to `na.rm` in [base::sum()], a logical value indicating
#'   whether NA values should be stripped before the computation proceeds.
#' @param align Character specifying whether the index of the result should be
#'   left- or right-aligned or centered (default) compared to the rolling window
#'   of observations.
#' @return The rolling sum.
#' @seealso [zoo::rollsum()]
#' @examples 
#' 
#' x <- 1:10
#' 
#' # Basic usage
#' rollsum(x, k = 3)
#' 
#' # Align at the right end, and use NA for padding
#' rollsum(x, k = 3, na_pad = TRUE, align = "right")
#' @export
rollsum <- function(x, k, na_pad = FALSE, na.rm = FALSE, align = c("center", "left", "right")) {
  align <- match.arg(align)
  
  align <- switch(align,
                  center = 1L,
                  left = 2L,
                  right = 3L)
  c_rollsum(x, k, na_pad = na_pad, na_rm = na.rm, align = align)
}


#' Calculate rolling means
#' 
#' Compared with [zoo::rollmean()], this function is implemented in C and is one
#' magnitude faster.
#' @param x A numeric vector.
#' @param k integer width of the rolling window.
#' @param na_pad Logical value about whether to pad left/right ends with NAs.
#' @param na.rm Similar to `na.rm` in [base::mean()], a logical value indicating
#'   whether NA values should be stripped before the computation proceeds.
#' @param align Character specifying whether the index of the result should be
#'   left- or right-aligned or centered (default) compared to the rolling window
#'   of observations.
#' @return The rolling mean.
#' @seealso [zoo::rollmean()]
#' @examples 
#' x <- 1:10
#' 
#' # Basic usage
#' rollmean(x, k = 3)
#' 
#' # Align at the right end, and use NA for padding
#' rollmean(x, k = 3, na_pad = TRUE, align = "right")
#' @export
rollmean <- function(x, k, na_pad = FALSE, na.rm = FALSE, align = c("center", "left", "right")) {
  align <- match.arg(align)
  
  align <- switch(align,
                  center = 1L,
                  left = 2L,
                  right = 3L)
  c_rollmean(x, k, na_pad = na_pad, na_rm = na.rm, align = align)
}