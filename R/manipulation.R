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
#'   functions applied on the corresponding columns, with the subset data frame
#'   as input ([data.table::.SD]). For example: `list(score = function(x)
#'   sum(x$score))` means for all intervals merged into one group, take the sum
#'   of the `score` column. Similar to `bedtools merge`'s `-c` and `-o`
#'   arguments.
#' @return A `data.table` object containing merged intervals.
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' merged <- merge_bed(bedtbl)
#' merged <- merge_bed(bedtbl, max_dist = 10, 
#'           operation = list(score1 = function(x) mean(x$score), 
#'                            score2 = function(x) sum(x$score)))
#' @references Manual page of `bedtools merge`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/merge.html}
#' @export
merge_bed <- function(x,
                  max_dist = 0,
                  operation = NULL) {
  stopifnot(max_dist >= 0)
  
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
    
    if (!is.null(operation)) {
      results2 <- names(operation) %>% map(function(op_name) {
        func <- operation[[op_name]]
        func(.SD)
      })
      names(results2) <- names(operation)
      c(results1, results2)
    } else
      results1
  }, by = c("chrom", (idx_colname))][, (idx_colname) := NULL]
  post_process_table(result)
}


# For each interval pair, get the overlapping start and end
overlapping_intervals <- function(x, y, overlap_table) {
  overlap_s <-
    as.integer(pmax(x[overlap_table$xid]$start - 1, y[overlap_table$yid]$start - 1))
  overlap_e <-
    as.integer(pmin(x[overlap_table$xid]$end, y[overlap_table$yid]$end))
  overlap_table[, c(.SD, list(overlap_s = overlap_s, overlap_e = overlap_e))]
}


# Filter `overlap_table` by `min_overlap`
# Return: a new `overlap_table` containing only the records pass the filter
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
#' @param x A `data.table` object.
#' @param y A `data.table` object.
#' @param mode Mode of the intersect operation. Can be one of the following:
#' 
#'   `default`: Similar to `bedtools intersect`'s default behavior.
#'
#'   `exclude`: Return features in `x` that do not overlap with `y`. Similar to
#'   `bedtools intersect`'s `-u` argument.
#'
#'   `wa`: Write the original entry in `x` for each overlap. Similar to
#'   `bedtools intersect`'s `-wa` argument.
#'
#'   `wb`: Write the original entry in `y` for each overlap. Useful for knowing
#'   what A overlaps. Similar to `bedtools intersect`'s `-wb` argument.
#'
#'   `wo`: Write the original `x` and `y` entries plus the number of base pairs
#'   of overlap between the two features. Only A features with overlap are
#'   reported.
#'
#'   `unique`: Write original `x` entry once if any overlaps found in `y`. In
#'   other words, just report the fact at least one overlap was found in `y`.
#'   Similar to `bedtools intersect`'s `-u` argument.
#'   
#'   `loj`: Perform a “left outer join”. That is, for each feature in `x` report
#'   each overlap with `j`. If no overlaps are found, report a NULL feature for B.
#' @param max_gap The largest gap for two intervals to be considered as
#'   overlapping. Default is 0 (no gap allowed).
#' @param min_overlap The smallest overlapping region for two intervals to be
#'   considered as overlapping. Default is 1.
#' @param min_overlap_type A character value indicating how `min_overlap` is
#'   interpreted. `bp` means `min_overlap` is the number of base pairs. `frac1`
#'   means `min_overlap` is the minimum overlap required as a fraction of `x`.
#'   Similarly, `frac2` means `min_overlap` is the minimum overlap required as a
#'   fraction of `y`. Similar to `bedtools intersect`'s `-f` and `-F` arguments.
#' @return A `data.table` representing the intersection of `x` and `y`
#' @examples 
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"))
#' 
#' # Basic usages
#' result <- intersect_bed(tbl_x, tbl_y)
#' 
#' # Exclude regions defined by tbl_y from tbl_x
#' result <- intersect_bed(tbl_x, tbl_y, mode = "exclude")
#' 
#' # For each overlap, return the original entries in tbl_x. For a interval in
#' # tbl_x, it is considered as overlapping only if 50% of it overlaps with an
#' # interval in tbl_y.
#' result <- intersect_bed(tbl_x, tbl_y, mode = "wa", min_overlap = 0.5, min_overlap_type = "frac1")
#' 
#' # For each overlap, return the original entries in both tbl_x and tbl_y, plus
#' # the number of overlapping base pairs. The minimum range for two intervals to
#' # be considered as overlapping is 5bp
#' result <- intersect_bed(tbl_x, tbl_y, mode = "wa", min_overlap = 5, min_overlap_type = "bp")
#' @references Manual page of `bedtools intersect`: \url{https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}
#' @export
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


# Load chrom sizes from hg19 or hg38, or any user-provided data file
load_chrom_sizes <- function(ref_genome) {
  chrom_sizes_file <- switch(
    ref_genome,
    hg19 = system.file("extdata", "human_g1k_v37.chrom.sizes", package = "bedtorch"),
    hg38 = system.file("extdata", "hg38.chrom.sizes", package = "bedtorch"),
    ref_genome
  )
  result <- fread(
    chrom_sizes_file,
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "size"),
    colClasses = c("character", "integer")
  )
  result[, chrom := factor(chrom, levels = str_sort(chrom, numeric = TRUE))]
  setkey(result, "chrom")
  result[]
}


#' Slop intervals
#' 
#' Similar to `bedtools slop`, this function increases the size of each feature
#' in a feature file by a user-defined number of bases. Will restrict the
#' resizing to the size of the chromosome.
#' @param x A `data.table` object.
#' @param chrom_sizes Provide the chromosome sizes data. If `hg19` or `hg38`,
#'   the embedded chromosome sizes for hg19 or hg38 will be used. Otherwise, it
#'   need to be the path of a user-provided file. If `NULL`, the resizing will
#'   not be restricted to the size of the chromosome.
#' @param slop The number of base pairs to slop, or the percentage of the
#'   interval's width (refer to `bedtools slop`'s `-pct` argument).
#' @param slop_type If `bp`, `slop` is an integer indicating the number of base
#'   pairs. If `frac`, `slop` is a float point number indicating the percentage.
#' @param mode If `both`, the resizing takes place on both ends of the interval.
#'   Otherwise, the resizing is only for the `left` end or the `right` end.
#' @return A slopped `data.tabe` object.
#' @examples 
#' # Load data
#' tbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' 
#' # Slop by 10 on both ends
#' result <- slop_bed(tbl, slop = 10L, slop_type = "bp", mode = "both")
#' 
#' # Slop by 10 only on the right end
#' result <- slop_bed(tbl, slop = 10L, slop_type = "bp", mode = "right")
#' 
#' # Slop by 20% on the left end
#' result <- slop_bed(tbl, slop = 0.2, slop_type = "frac", mode = "left")
#' @references Manual page of `bedtools slop`: \url{https://bedtools.readthedocs.io/en/latest/content/tools/slop.html}
#' @export
slop_bed <-
  function(x,
           chrom_sizes = c("hg19", "hg38"),
           slop = 1L,
           slop_type = c("bp", "frac"),
           mode = c("both", "left", "right")) {
    mode <- match.arg(mode)
    slop_type <- match.arg(slop_type)
    chrom_sizes <- if (is.null(chrom_sizes))
      NULL
    else {
      chrom_sizes <- tryCatch({
        match.arg(chrom_sizes)
      }, error = function(e) {return(chrom_sizes)})
      load_chrom_sizes(chrom_sizes)
    }
      
    x <- data.table::copy(x)
    
    # Make sure each interval has a positive width after operation
    if (slop_type == "frac") {
      stopifnot((mode == "both" &&
                   slop > -0.5) || (mode != "both" && slop > -1))
    } else if (slop < 0) {
      if (mode == "both")
        negative_idx <- which(x[, end - start] < 2 * abs(slop))
      else
        negative_idx <- which(x[, end - start] < abs(slop))
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
    
    if (slop_type == "frac") {
      width_name <- .available_colname(x, "width")
      x[, (width_name) := end - start]
      on.exit({
        x[, (width_name) := NULL]
      }, add = TRUE)
      
      if (mode == "both")
        x[, `:=`(start = start - floor(slop * width),
                 end = end + floor(slop * width))]
      else if (mode == "left")
        x[, start := start - floor(slop * width)]
      else if (mode == "right")
        x[, end := end + floor(slop * width)]
      else
        stop(paste0("Invalid mode: ", mode))
    } else if (slop_type == "bp") {
      if (mode == "both")
        x[, `:=`(start = start - slop, end = end + slop)]
      else if (mode == "left")
        x[, start := start - slop]
      else if (mode == "right")
        x[, end := end + slop]
      else
        stop(paste0("Invalid mode: ", mode))
    } else
      stop(paste0("Invalid slop type: ", slop_type))
    
    x[start < 0, start := 0]
    
    if (!is.null(chrom_sizes)) {
      n1 <- nrow(x)
      idx <- x[, 1:3][chrom_sizes, nomatch=0][, id := 1:.N]
      n2 <- nrow(idx)
      
      if (n2 < n1)
        warning(str_interp("${n1 - n2} out of ${n1} rows: cannot find corresponding chroms in the reference genome"))
      
      # Where the slopped interval exceeds the chrom size limit
      idx <- idx[end > size]
      x[idx$id, end := idx$size]
    }
    
    x[, `:=`(start = as.integer(start), end = as.integer(end))]
    setkey(x, "chrom", "start", "end")
    x[]
  }


#' Map over scaffold intervals
#'
#' Map overlapping features in `data` onto intervals in `scaffold` and apply
#' statistics and/or summary operations on those features.
#' @param data A `data.table` object.
#' @param scaffold A `data.table` object containing intervals upon which you
#'   want to map `data`.
#' @param operation List of functions for the statistics and summary operations.
#'   This is similar to [bedtorch::merge_bed()]
#' @param min_overlap See [bedtorch::intersect_bed()].
#' @param min_overlap_type [bedtorch::intersect_bed()].
#' @return A mapped `data.table` object.
#' @examples
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"))
#'
#' # Basic usage
#' result <- map_bed(tbl_x, tbl_y, operation = list(score_mean = function(x) mean(x$score)))
#'
#' # Perform the mapping, requiring the minimum overlapping of 5bp
#' result <- map_bed(tbl_x, tbl_y, operation = list(score_mean = function(x) mean(x$score)), 
#'                   min_overlap = 5, min_overlap_type = "bp")
#' @references Manual page of `bedtools map`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#' @seealso [bedtorch::merge_bed()], [bedtorch::intersect_bed()]
#' @export
map_bed <- function(data, scaffold, operation, 
                    min_overlap = 1,
                    min_overlap_type = c("bp", "frac1", "frac2")) {
  stopifnot(!is.null(operation))
  
  result <-
    intersect_bed(
      data,
      scaffold[, 1:3],
      mode = "wo",
      min_overlap = min_overlap,
      min_overlap_type = min_overlap_type
    )
  
  # chrom/start/end from scaffold now becomes the main chrom/start/end
  fields <- colnames(result) %>% tail(-3) %>% head(-1)
  fields <- c(tail(fields, 3), head(fields, -3))
  result <- result[, ..fields]
  setnames(result,
           old = 1:3, 
           new = c("chrom", "start", "end"))
  setkey(result, "chrom", "start", "end")
  
  result[, {
    # Apply functions
    dt2 <- names(operation) %>% map(function(op_name) {
      func <- operation[[op_name]]
      func(.SD)
    })
    names(dt2) <- names(operation)
    dt2
  }, by = c("chrom", "start", "end")][scaffold[, 1:3]]
}
