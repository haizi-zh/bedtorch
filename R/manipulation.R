#' Merge overlapping intervals
#' 
#' This operation is similar to `bedtools merge`.
#' @param x A `data.table` object.
#' @param max_dist Maximum distance between features allowed for features to be
#'   merged. Default is 0. That is, overlapping and/or book-ended features are
#'   merged.
#' @param operation Specify what operations should be applied to merged
#'   intervals. Default is `NULL`, i.e. do not apply any operation and only
#'   return the first three columns (`chrom`, `start`, `end`). Must be a list in
#'   the following format: `list(col1 = func1, col2 = func2, ...)`, where `col1`
#'   and `col2` are column names, and `func1` and `func2` are univariate
#'   functions applied on the corresponding columns. For example: `list(score =
#'   sum)` means for all intervals merged into one group, take the sum of the
#'   `score` column. Similar to `bedtools merge`'s `-c` and `-o` arguments.
#' @examples merge_bed(x, operation = list(score = mean))
#' @seealso https://bedtools.readthedocs.io/en/latest/content/tools/merge.html
#' @return A `data.table` object containing merged intervals.
merge_bed <- function(x,
                  max_dist = 0,
                  operation = NULL) {
  stopifnot(max_dist >= 0)
  if (!is.null(operation))
    stopifnot(all(names(operation) %in% colnames(x)))
  
  idx_colname <- .available_colname(x, "idx")
  x[, (idx_colname) := head(c(0, cumsum(shift(start, n = -1) > cummax(end) + max_dist)), n = -1), by = chrom]
  on.exit({
    x[, (idx_colname) := NULL]
  }, add = TRUE)
  
  result <- x[, {
    start = min(start)
    end = max(end)
    
    results1 <- list(start = start, end = end)
    # Apply functions
    r <- {
      if (!is.null(operation)) {
        results2 <- names(operation) %>% map(function(field) {
          func <- operation[[field]]
          func(.SD[[field]])
        })
        names(results2) <- names(operation)
        c(results1, results2)
      } else
        results1
    }
    r
  }, by = c("chrom", (idx_colname))][, (idx_colname) := NULL]
  post_process_table(result)
}


#' For each interval pair, get the overlapping start and end
#' 
overlapping_intervals <- function(x, y, overlap_table) {
  overlap_s <-
    as.integer(pmax(x[overlap_table$xid]$start - 1, y[overlap_table$yid]$start - 1))
  overlap_e <-
    as.integer(pmin(x[overlap_table$xid]$end, y[overlap_table$yid]$end))
  overlap_table[, c(.SD, list(overlap_s = overlap_s, overlap_e = overlap_e))]
}


#' Filter `overlap_table` by `min_overlap`
#' 
#' @return A new `overlap_table` containing only the records pass the filter
min_overlap_filter <- function(x, y, overlap_table, min_overlap, min_overlap_type) {
  if (min_overlap_type == "bp" && min_overlap > 1) {
    overlap_table <-
      overlapping_intervals(x, y, overlap_table = overlap_table)
    overlap_table <-
      overlap_table[overlap_e - overlap_s >= min_overlap]
  } else if (min_overlap_type == "frac1") {
    overlap_table <-
      overlapping_intervals(x, y, overlap_table = overlap_table)
    interval_width <- x[overlap_table$xid][, end - start + 1]
    overlap_table <-
      overlap_table[overlap_e - overlap_s >= ceiling(interval_width * min_overlap)]
  } else if (min_overlap_type == "frac2") {
    overlap_table <-
      overlapping_intervals(x, y, overlap_table = overlap_table)
    interval_width <- y[overlap_table$yid][, end - start + 1]
    overlap_table <-
      overlap_table[overlap_e - overlap_s >= ceiling(interval_width * min_overlap)]
  } else {
    stop(str_intert("Invalid minimal overlap type: ${min_overlap_type}"))
  }
  overlap_table
}


#' Apply `intersect` operation between two tables
#'
#' @param mode `default`: Similar to `bedtools intersect`'s default behavior.
#'
#'   `wa`: Similar to `bedtools intersect`'s `-wa` argument.
#'
#'   `wb`: Similar to `bedtools intersect`'s `-wb` argument.
#'
#'   `unique`: Similar to `bedtools intersect`'s `-u` argument.
#' @param max_gap The largest gap for two intervals to be considered as
#'   overlapping. Default is 0 (no gap allowed).
#' @param min_overalp The smallest overlapping region for two intervals to be
#'   considered as overlapping. Default is 1.
intersect_bed <-
  function(x,
           y,
           mode = c("default", "exclude", "wa", "wb", "wo", "unique", "loj"),
           max_gap = 0,
           min_overlap = 1,
           min_overlap_type = c("bp", "frac1", "frac2")) {
  mode <- match.arg(mode)
  min_overlap_type <- match.arg(min_overlap_type)
  
  stopifnot(max_gap >= 0)
  stopifnot(min_overlap > 0)
  
  if (max_gap > 0) stop("max_gap has not been implemented yet.")
  
  # Change it to 1-based convention
  x[, start := start + 1L]
  y[, start := start + 1L]
  setkey(x, "chrom", "start", "end")
  setkey(y, "chrom", "start", "end")
  on.exit({
    x[, start := start - 1L]
    y[, start := start - 1L]
  }, add = TRUE)
  
  if (max_gap != 0) {
    y[, `:=`(start = start - max_gap, end = end + max_gap)]
    setkey(y, "chrom", "start", "end")
    on.exit({
      y[, `:=`(start = start + max_gap, end = end - max_gap)]
    }, add = TRUE)
  }
  
  on.exit({
    setkey(x, "chrom", "start", "end")
    setkey(y, "chrom", "start", "end")
  }, add = TRUE)
  
  if (mode == "exclude") {
    result <- data.table::foverlaps(
      x,
      y[, 1:3],
      by.x = c("chrom", "start", "end"),
      by.y = c("chrom", "start", "end"),
      nomatch = NA
    )
    
    exclude_idx <- result[, {
      overlap_s <- ifelse(i.start < start, start, i.start)
      overlap_e <- ifelse(i.end > end, end, i.end)
      is.na(start) | (overlap_e - overlap_s + 1 < min_overlap)
    }]
    result <- result[exclude_idx]
    result[, `:=`(start = NULL, end = NULL)]
    setnames(result, c("i.start", "i.end"), c("start", "end"))
    
    # Change back to 0-based
    result[, start := start - 1L]
    result[, `:=`(start = as.integer(start), end = as.integer(end))]
    setkey(result, "chrom", "start", "end")
    
    return(result)
  }
  
  # Mode: default, wa, etc.
  overlap_table <- data.table::foverlaps(
    x[, 1:3],
    y[, 1:3],
    by.x = c("chrom", "start", "end"),
    by.y = c("chrom", "start", "end"),
    nomatch = NULL,
    which = TRUE
  )
  
  # When we need to calculate the overlapping intervals
  # Some mode requires overlap interval data
  if ((min_overlap_type=="bp" && min_overlap > 1) ||
      (min_overlap_type %in% c("frac1", "frac2"))) {
    overlap_table <- min_overlap_filter(
      x = x,
      y = y,
      overlap_table = overlap_table,
      min_overlap = min_overlap,
      min_overlap_type = min_overlap_type
    )
  } else if (mode %in% c("default", "wo", "wb")) {
    # No need to run the min_overlap filter, but need the overlap interval data
    overlap_table <- overlapping_intervals(x, y, overlap_table = overlap_table)
  }
  
  if (mode %in% c("default", "wb")) {
    result <- x[overlap_table$xid]
    # Change back to 0-based
    result[, `:=`(start = overlap_table$overlap_s, end = overlap_table$overlap_e)]
    
    if (mode == "wb") {
      # Column binding
      result2 <- y[overlap_table$yid]
      result2[, start := start - 1L]
      
      result <- setDT(
        unlist(list(result, result2), recursive = FALSE),
        check.names = TRUE
      )
      # # # of overlapping base pairs
      # overlap_col <- .available_colname(result, "overlap")
      # result[, (overlap_col) := as.integer(overlap_table$overlap_e - overlap_table$overlap_s)]
    }
    
    setkey(result, "chrom", "start", "end")
    return(result)
  }
  
  if (mode %in% c("unique", "wa", "wo")) {
    result <- x[overlap_table$xid]
    # Change back to 0-based
    result[, start := start - 1L]
    if (mode == "wo") {
      # Column binding
      result2 <- y[overlap_table$yid]
      result2[, start := start - 1L]
      
      result <- setDT(
        unlist(list(result, result2), recursive = FALSE),
        check.names = TRUE
      )
      # # of overlapping base pairs
      overlap_col <- .available_colname(result, "overlap")
      result[, (overlap_col) := as.integer(overlap_table$overlap_e - overlap_table$overlap_s)]
    } else if (mode == "unique") {
      result <- unique(result)
    }
    setkey(result, "chrom", "start", "end")
    return(result)
  }
  
  stop(str_interp("${mode} not implemented yet"))
}



#' Grow intervals
#' 
#' This function resizes each interval. 
grow_bed <- function(x, grow = 1L, grow_type = c("bp", "frac"), mode = c("both", "left", "right")) {
  mode <- match.arg(mode)
  
  # Make sure each interval has a positive width after operation
  if (grow_type == "frac") {
    stopifnot((mode == "both" &&
                 grow > -0.5) || (mode != "both" && grow > -1))
  } else if (grow < 0) {
    if (mode == "both")
      negative_idx <- which(x[, end - start] < 2 * abs(grow))
    else
      negative_idx <- which(x[, end - start] < abs(grow))
    if (length(negative_idx) > 0) {
      # Some intervals are ne
      interval <- x[negative_idx[1]]
      stop(
        str_interp(
          "Interval ${interval$chrom}:${interval$start + 1}-${interval$end} becomes negative after operation"
        )
      )
    }
  }
  
  if (grow_type == "frac") {
    width_name <- .available_colname(x, "width")
    x[, (width_name) := end - start]
    on.exit({
      x[, (width_name) := NULL]
    }, add = TRUE)
    
    if (mode == "both") 
      x[, `:=`(start = round(start - grow * width), end = round(end + grow * width))]
    else if (mode == "left")
      x[, start := round(start - grow * width)]
    else if (mode == "right")
      x[, end := round(end + grow * width)]
    else
      stop(paste0("Invalid mode: ", mode))
  } else if (grow_type == "bp") {
    if (mode == "both")
      x[, `:=`(start = start - grow, end = end + grow)]
    else if (mode == "left")
      x[, start := start - grow]
    else if (mode == "right")
      x[, end := end + grow]
    else
      stop(paste0("Invalid mode: ", mode))
  } else
      stop(paste0("Invalid grow type: ", grow_type))
  
  x[start < 0, start := 0]
    
  x[, `:=`(start = as.integer(start), end = as.integer(end))]
  setkey(x, "chrom", "start", "end")
  x[]
}