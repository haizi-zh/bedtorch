# Find a new column name which does not conflict with existing names
# 
# # Returns "chrom.1"
# .available_colname(dt, "chrom")
.available_colname <- function(dt, name) {
  if (is(dt, "GRanges"))
    fields <- colnames(mcols(dt))
  else
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


# Constructor of bedtorch_table from data.frame
# 
# @param dt A `data.table` or `data.frame` input.
# @param genome A character value specifying the genome name.
new_bedtorch_table <- function(dt, genome = NULL) {
  if (is.null(dt))
    return(dt)

  if (is(dt, "bedtorch_table"))
    return(dt)
  
  stopifnot(is.null(genome) || (is.character(genome) && length(genome) == 1))
  
  if (is(dt, "data.frame") && !is(dt, "data.table"))
    data.table::setDT(dt)
  
  # data.table requires modify-in-place
  dt_classes <- class(dt)
  # Strip leading "bedtorch_table" entries, in case there are not duplicates
  non_matches <- which(dt_classes != "bedtorch_table")
  # rare case, where dt_classes all repeative "bedtorch_table"
  # DO NOT TOUCH
  stopifnot(length(non_matches) > 0)
  # In case dt is alread bedtorch_table
  dt_classes <- c("bedtorch_table", dt_classes[min(non_matches):length(dt_classes)])
  data.table::setattr(dt, "class", dt_classes)
  
  # Check first three columns
  cols <- colnames(dt)
  if (length(cols) == 0)
    return(dt)
  
  stopifnot(length(cols) >= 3)
  data.table::setnames(dt, 1:3, c("chrom", "start", "end"))
  
  data.table::setkey(dt, "chrom", "start", "end")
  if (!is.null(genome))
    data.table::setattr(dt, "genome", genome)
  else
    data.table::setattr(dt, "genome", NULL)
  dt
}


#' Convert a `data.frame` representation to `GenomicRanges` representation
#' 
#' @param x An input `data.frame`.
#' @return `GenomicRanges` converted from `x`.
#' @examples
#' dt <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), 
#'                use_gr = FALSE, genome = "hs37-1kg")
#' as.GenomicRanges(dt)
#' @export
as.GenomicRanges <- function(x) {
  UseMethod("as.GenomicRanges")
}

#' @export
as.GenomicRanges.GRanges <- function(x) {
  return(x)
}


#' @export
as.GenomicRanges.data.frame <- function(x) {
  dt <- x

  genome <- attr(dt, "genome")
  stopifnot(is.null(genome) || (is.character(genome) && length(genome) == 1))
  
  if (is.null(x))
    NULL
  else if (nrow(x) == 0)
    GenomicRanges::GRanges()
  else
    GenomicRanges::makeGRangesFromDataFrame(
      dt,
      keep.extra.columns = TRUE,
      seqinfo = get_seqinfo(genome),
      starts.in.df.are.0based = TRUE
    )
}


#' Convert a `GenomicRanges` or a `data.frame` BED to standard `data.table` BED
#' 
#' @param x An input. Must be either `GenomicRanges` or `data.frame`. If a
#'   `data.frame`, the first three columns should be chrom, start and end. For
#'   column names, refer to [GenomicRanges::makeGRangesFromDataFrame()].
#' @return `data.table` converted from `x`.
#' @export
#' @examples
#' gr <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), 
#'                use_gr = TRUE, genome = "hs37-1kg")
#' as.bedtorch_table(gr)
as.bedtorch_table <- function(x, genome = NULL) {
  UseMethod("as.bedtorch_table")
}


#' @export
as.bedtorch_table.bedtorch_table <- function(x, genome = NULL) {
  if (is.null(genome))
    return(x)
  else {
    data.table::setattr(x, "genome", genome)
    return(x)
  }
}


#' @export
as.bedtorch_table.data.table <- function(x, genome = NULL) {
  new_bedtorch_table(x, genome = genome %||% attr(x, "genome"))
}


#' @export
as.bedtorch_table.data.frame <- function(x, genome = NULL) {
  data.table::setDT(x)
  new_bedtorch_table(x, genome = genome %||% attr(x, "genome"))
}


#' @export
as.bedtorch_table.GRanges <- function(x) {
  gr <- x
  dt <- data.table::as.data.table(gr)
  data.table::setnames(dt, "seqnames", "chrom")
  dt[, `:=`(start = as.integer(start - 1), width = NULL, strand = NULL)]
  data.table::setkey(dt, "chrom", "start", "end")
  
  genome <- GenomeInfoDb::genome(gr) %>% unique()
  stopifnot(length(genome) == 1)
  if (is.na(genome))
    genome <- NULL
  
  new_bedtorch_table(dt, genome = genome)
}


#' @export
print.bedtorch_table <- function(x, ...) {
  NextMethod()
  
  genome <- attr(x, "genome")
  cat("-------\n")
  if (is.null(genome))
    genome <- "unspecified"
  cat(str_interp("genome: ${genome}.\n"))
}

