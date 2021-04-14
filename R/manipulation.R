#' Merge overlapping intervals
#' 
#' This operation is similar to `bedtools merge`.
#' @param x A `GRanges` object.
#' @param max_dist Maximum distance between features allowed for features to be
#'   merged. Default is 0. That is, overlapping and/or book-ended features are
#'   merged.
#' @param operation Specify what operations should be applied to merged
#'   intervals. Default is `NULL`, i.e. do not apply any operation and only
#'   return the first three columns (`chrom`, `start`, `end`). Must be a list in
#'   the following format: 
#'   `list(col1_save = list(on = col1, func = func1), col2_save = list(on = col2, func = func2), ...)`, 
#'   where `col1` and `col2` are column names to apply merge functions,
#'   `col1_save` and `col2_save` are columns to save the results, and `func1`
#'   and `func2` are univariate merge functions. For example: `list(max_score =
#'   list(on = "score", func = max))` means for all intervals merged into one
#'   group, take the sum of the `score` column, and save to `max_score`. Similar
#'   to `bedtools merge`'s `-c` and `-o` arguments.
#' @return A `GRanges` object containing merged intervals.
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' merged <- merge_bed(bedtbl)
#' head(merged)
#' 
#' merged <- merge_bed(bedtbl, max_dist = 10, 
#'           operation = list(score1 = list(on = "score", func = mean),
#'                            score2 = list(on = "score", func = sum)))
#' head(merged)
#' @references Manual page of `bedtools merge`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/merge.html}
#' @export
merge_bed <- function(x,
                      max_dist = 0,
                      operation = NULL) {
  UseMethod("merge_bed")
}

#' @export
merge_bed.bedtorch_table <- function(x,
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
      results2 <- lapply(names(operation), function(op_name) {
        func <- operation[[op_name]]
        func(.SD[[op_name]])
      })
      names(results2) <- names(operation)
      c(results1, results2)
    } else
      results1
  }, by = c("chrom", (idx_colname))][, (idx_colname) := NULL]
  post_process_table(result)
  
  new_bedtorch_table(result, genome = attr(x, "genome"))
}

#' @export
merge_bed.GRanges <- function(x, max_dist = 0, operation = NULL) {
  stopifnot(max_dist >= 0)
  
  merged <- GenomicRanges::reduce(x,
                                  min.gapwidth = max_dist + 1L,
                                  with.revmap = !is.null(operation))
  
  if (!is.null(operation)) {
    aggregated <- names(operation) %>% map(function(op_name) {
      func <- operation[[op_name]]$func
      col_name <- operation[[op_name]]$on
      data <- mcols(x)[[col_name]]
      # func(IRanges::extractList(data, merged$revmap))
      sapply(IRanges::extractList(data, merged$revmap), func)
    })
    names(aggregated) <- names(operation)
    mcols(merged) <- c(mcols(merged), aggregated)
  }
  
  if (!is.null(operation))
    merged$revmap = NULL
  
  merged
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
#' @param x A `GRanges` object.
#' @param y A `GRanges` object.
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
#'   overlapping. Default is -1 (no gap allowed, adjacent intervals not allowd).
#' @param min_overlap The smallest overlapping region for two intervals to be
#'   considered as overlapping. Default is 0.
#' @param min_overlap_type A character value indicating how `min_overlap` is
#'   interpreted. `bp` means `min_overlap` is the number of base pairs. `frac1`
#'   means `min_overlap` is the minimum overlap required as a fraction of `x`.
#'   Similarly, `frac2` means `min_overlap` is the minimum overlap required as a
#'   fraction of `y`. Similar to `bedtools intersect`'s `-f` and `-F` arguments.
#' @return A `GRanges` representing the intersection of `x` and `y`
#' @examples 
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"), use_gr = FALSE)
#' 
#' # Basic usages
#' result <- intersect_bed(tbl_x, tbl_y)
#' head(result)
#' 
#' # Exclude regions defined by tbl_y from tbl_x
#' result <- intersect_bed(tbl_x, tbl_y, mode = "exclude")
#' head(result)
#' 
#' # For each overlap, return the original entries in tbl_x. For a interval in
#' # tbl_x, it is considered as overlapping only if 50% of it overlaps with an
#' # interval in tbl_y.
#' result <- intersect_bed(tbl_x, tbl_y, mode = "wa", min_overlap = 0.5, min_overlap_type = "frac1")
#' head(result)
#' 
#' # For each overlap, return the original entries in both tbl_x and tbl_y, plus
#' # the number of overlapping base pairs. The minimum range for two intervals to
#' # be considered as overlapping is 5bp
#' result <- intersect_bed(tbl_x, tbl_y, mode = "wa", min_overlap = 5, min_overlap_type = "bp")
#' head(result)
#' @references Manual page of `bedtools intersect`: \url{https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}
#' @seealso [GenomicRanges::findOverlaps()]
#' @export
intersect_bed <-
  function(x,
           y,
           mode = c("default", "exclude", "wa", "wb", "wo", "unique", "loj"),
           max_gap = -1L,
           min_overlap = 0L,
           min_overlap_type = c("bp", "frac1", "frac2")) {
    UseMethod("intersect_bed")
  }


#' @export
intersect_bed.bedtorch_table <-
  function(x,
           y,
           mode = c("default", "exclude", "wa", "wb", "wo", "unique", "loj"),
           max_gap = -1L,
           min_overlap = 0L,
           min_overlap_type = c("bp", "frac1", "frac2")) {
    stopifnot(is(x, "bedtorch_table"))
    stopifnot(is(y, "bedtorch_table"))
  mode <- match.arg(mode)
  min_overlap_type <- match.arg(min_overlap_type)
  
  stopifnot(max_gap >= -1L)
  stopifnot(min_overlap >= 0L)
  
  # Determine the genome attribute of the reuslt
  genome <- local({
    genome_list <-
      list(x, y) %>% map( ~ attr(., which = "genome")) %>% 
      discard(is.null) %>% unique()
    
    if (length(genome_list) == 2)
      # Dual genome, can't decide the genome of the result
      NULL
    else if (length(genome_list) == 0)
      NULL
    else
      genome_list[[1]]
  })
  
  if (max_gap != -1L) stop("max_gap has not been implemented yet.")
  
  # Here's a trick: for GenomicRanges::findOverlaps,
  # when 'type' is "any", at least one of 'maxgap' and 'minoverlap' must be set
  # to its default value
  # By default, max_gap = -1 means STRICT inside the other, which implies positive
  # overlap even if min_overlap = 0
  # This works for GenomicRanges, but not data.table, since for data.table max_gap
  # has not been implemented yet. Thus max_gap is simply ignored. As a result, we
  # shouldn't allow min_overlap to be 0
  if (min_overlap == 0)
    min_overlap = 1L
  
  # Change it to 1-based convention
  x[, start := start + 1L]
  y[, start := start + 1L]
  setkey(x, "chrom", "start", "end")
  setkey(y, "chrom", "start", "end")
  on.exit({
    x[, start := start - 1L]
    y[, start := start - 1L]
  }, add = TRUE)
  
  # if (max_gap != 0) {
  #   y[, `:=`(start = start - max_gap, end = end + max_gap)]
  #   setkey(y, "chrom", "start", "end")
  #   on.exit({
  #     y[, `:=`(start = start + max_gap, end = end - max_gap)]
  #   }, add = TRUE)
  # }
  
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
    
    return(new_bedtorch_table(result, genome))
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
    return(new_bedtorch_table(result, genome))
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
    return(new_bedtorch_table(result, genome))
  }
  
  stop(str_interp("${mode} not implemented yet"))
}


intersect_bed_unique <- function(x, y, max_gap = -1L, min_overlap = 0L) {
  IRanges::subsetByOverlaps(x, y, maxgap = max_gap, minoverlap = min_overlap)
}


intersect_bed_default <- function(x, y, max_gap = -1L, min_overlap = 0L) {
  stopifnot(max_gap == -1L)

  hits <- findOverlaps(x, y, minoverlap = min_overlap)
  
  result <- x[from(hits)]
  ranges(result) <- pintersect(ranges(x)[from(hits)], ranges(y)[to(hits)])
  result
}


#' @export
intersect_bed.GRanges <-
  function(x,
           y,
           mode = c("default", "exclude", "wa", "wb", "wo", "unique", "loj"),
           max_gap = -1L,
           min_overlap = 0L,
           min_overlap_type = c("bp", "frac1", "frac2")) {
    mode <- match.arg(mode)
    min_overlap_type <- match.arg(min_overlap_type)
    
    if (min_overlap_type != "bp")
      stop("Currently min_overlap_type other than bp is not supported")
    
    if (mode == "unique")
      intersect_bed_unique(x, y, max_gap = max_gap, min_overlap = min_overlap)
    else if (mode == "default") {
      if (max_gap != -1L)
        stop("In default mode, max_gap must be the default value (-1)")
      intersect_bed_default(x, y, max_gap = max_gap, min_overlap = min_overlap)
    } else
      stop(str_interp("Unsupported mode ${mode}"))
  }


#' Exclude certain regions
#' 
#' Return original entries in `x` that do not overlap with any intervals in `y`.
#' Essentially, this function is just a short-hand for
#' [bedtorch::intersect_bed()] in `exclude` mode.
#' @param x A `GRanges` object.
#' @param y A `GRanges` object.
#' @param max_gap The largest gap for two intervals to be considered as
#'   overlapping. Default is -1 (no gap allowed, adjacent intervals not allowd).
#' @param min_overlap The smallest overlapping region for two intervals to be
#'   considered as overlapping. Default is 0.
#' @param min_overlap_type A character value indicating how `min_overlap` is
#'   interpreted. `bp` means `min_overlap` is the number of base pairs. `frac1`
#'   means `min_overlap` is the minimum overlap required as a fraction of `x`.
#'   Similarly, `frac2` means `min_overlap` is the minimum overlap required as a
#'   fraction of `y`. Similar to `bedtools intersect`'s `-f` and `-F` arguments.
#' @return A `GRanges`.
#' @examples 
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"), use_gr = FALSE)
#' 
#' # Basic usages
#' head(exclude_bed(tbl_x, tbl_y))
#' @seealso [bedtorch::intersect_bed()] [GenomicRanges::findOverlaps()]
#' @references Manual page of `bedtools intersect`: \url{https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html}
#' @export
exclude_bed <- function(x,
                        y,
                        max_gap = -1L,
                        min_overlap = 0L,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  UseMethod("exclude_bed")
}


#' @export
exclude_bed.GRanges <- function(x,
                        y,
                        max_gap = -1L,
                        min_overlap = 0L,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  min_overlap_type <- match.arg(min_overlap_type)
  if (min_overlap_type != "bp")
    stop("Currently min_overlap_type other than bp is not supported")
  
  hits <- findOverlaps(x, y, maxgap = max_gap, minoverlap = min_overlap)
  
  # Non-hit intervals
  all_ids <- seq.int(length(x))
  hit_ids <- from(hits)
  
  x[sort(setdiff(all_ids, hit_ids))]
}


#' @export
exclude_bed.bedtorch_table <- function(x,
                        y,
                        max_gap = -1L,
                        min_overlap = 0L,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  stopifnot(is(x, "bedtorch_table"))
  stopifnot(is(y, "bedtorch_table"))
  intersect_bed.bedtorch_table(
    x,
    y,
    mode = "exclude",
    max_gap = max_gap,
    min_overlap = min_overlap,
    min_overlap_type = min_overlap_type
  )
}


# Load chrom sizes, or any user-provided data file
load_chrom_sizes <- function(genome) {
  seqinfo <- get_seqinfo(genome)
  chrom_levels <- str_sort(unique(seqnames(seqinfo)), numeric = TRUE)
  chrom_sizes <- data.table(
    chrom = factor(seqnames(seqinfo), levels = chrom_levels),
    size = seqlengths(seqinfo) %>% unname()
  )
  setkey(chrom_sizes, "chrom")
  chrom_sizes
}
  

#' Slop intervals
#' 
#' Similar to `bedtools slop`, this function increases the size of each feature
#' in a feature file by a user-defined number of bases. Will restrict the
#' resizing to the size of the chromosome.
#' @param x A `GRanges` object.
#' @param genome Specify the reference genome for the BED file. `genome` can be
#'   a valid genome name in [GenomeInfoDb::Seqinfo()], e.g. `GRCh37`, or
#'   `hs37-1kg`, which is a genome shipped with this package, or any custom
#'   chromosome size files (local or remote). If `NULL`, will try to obtain such
#'   information from input data. Refer to [bedtorch::read_bed()].
#' @param slop The number of base pairs to slop, or the percentage of the
#'   interval's width (refer to `bedtools slop`'s `-pct` argument).
#' @param slop_type If `bp`, `slop` is an integer indicating the number of base
#'   pairs. If `frac`, `slop` is a float point number indicating the percentage.
#' @param mode If `both`, the resizing takes place on both ends of the interval.
#'   Otherwise, the resizing is only for the `left` end or the `right` end.
#' @return A slopped `GRanges` object.
#' @examples 
#' # Load data
#' tbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), 
#'                 use_gr = FALSE, genome = "hs37-1kg")
#' 
#' # Slop by 10 on both ends
#' result <- slop_bed(tbl, slop = 10L, slop_type = "bp", mode = "both")
#' 
#' # Slop by 10 only on the right end
#' result <- slop_bed(tbl, slop = 10L, slop_type = "bp", mode = "right")
#' head(result)
#' 
#' # Slop by 20% on the left end
#' result <- slop_bed(tbl, slop = 0.2, slop_type = "frac", mode = "left")
#' head(result)
#' @references Manual page of `bedtools slop`: \url{https://bedtools.readthedocs.io/en/latest/content/tools/slop.html}
#' @export
slop_bed <-
  function(x,
           genome = NULL,
           slop = 1L,
           slop_type = c("bp", "frac"),
           mode = c("both", "left", "right")) {
    UseMethod("slop_bed")
  }


#' @export
slop_bed.GRanges <-
  function(x,
           genome = NULL,
           slop = 1L,
           slop_type = c("bp", "frac"),
           mode = c("both", "left", "right")) {
    dt <- as.bedtorch_table(x)
    slop_bed(
      dt,
      genome = genome,
      slop = slop,
      slop_type = slop_type,
      mode = mode
    ) %>%
      as.GenomicRanges()
  }


#' @export
slop_bed.bedtorch_table <-
  function(x,
           genome = NULL,
           slop = 1L,
           slop_type = c("bp", "frac"),
           mode = c("both", "left", "right")) {
    mode <- match.arg(mode)
    slop_type <- match.arg(slop_type)
    
    if (is.null(genome))
      genome <- attr(x, "genome")
    if (is.null(genome))
      stop("Must provide genome information")
    
    chrom_sizes <- load_chrom_sizes(genome)
      
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
      idx <- x[, 1:3][chrom_sizes, nomatch=0][, id := seq_len(.N)]
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

#' @export
map_bed.GRanges <- function(x, y, operation, 
                    max_gap = -1L,
                    min_overlap = 0L,
                    min_overlap_type = c("bp", "frac1", "frac2")) {
  min_overlap_type <- match.arg(min_overlap_type)
  
  stopifnot(all(operation %>% map_chr(~ .$on) %in% names(mcols(x))))
  
  if (min_overlap_type != "bp")
    stop("Currently min_overlap_type other than bp is not supported")
  
  hits <- findOverlaps(x, y, maxgap = max_gap, minoverlap = min_overlap)
  
  # Use data.table as a helper
  col_names <- operation %>% map_chr(~ .$on)
  dt <- as.data.table(mcols(x[from(hits)])[col_names])[, hits_id := to(hits)]

  aggregated <- names(operation) %>% map(function(op_name) {
    func <- operation[[op_name]]$func
    col_name <- operation[[op_name]]$on
    dt[, func(.SD[[col_name]]), by = hits_id]$V1
  })
  names(aggregated) <- names(operation)

  result <- y[unique(to(hits))]
  mcols(result) <- aggregated
  result
  
  # # Stick to GenomicRanges
  # result <- y[unique(to(hits))]
  # 
  # aggregated <- names(operation) %>% map(function(op_name) {
  #   func <- operation[[op_name]]$func
  #   col_name <- operation[[op_name]]$on
  #   v <- aggregate(mcols(x[from(hits)])[col_name], list(idx = to(hits)), func)
  #   v[[col_name]]
  # })
  # names(aggregated) <- names(operation)
  # 
  # mcols(result) <- aggregated
  # result
}


#' Map over scaffold intervals
#'
#' Map overlapping features in `data` onto intervals in `scaffold` and apply
#' statistics and/or summary operations on those features.
#' @param x A `GRanges` object.
#' @param y A `GRanges` object containing intervals upon which you
#'   want to map `data`.
#' @param operation List of functions for the statistics and summary operations.
#'   This is similar to [bedtorch::merge_bed()]
#' @param max_gap The largest gap for two intervals to be considered as
#'   overlapping. Default is -1 (no gap allowed, adjacent intervals not allowd).
#' @param min_overlap See [bedtorch::intersect_bed()].
#' @param min_overlap_type [bedtorch::intersect_bed()].
#' @return A mapped `GRanges` object.
#' @examples
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"))
#'
#' # Basic usage
#' result <- map_bed(tbl_x, tbl_y, operation = list(score_mean = list(on = "score", func = mean)))
#' head(result)
#'
#' # Perform the mapping, requiring the minimum overlapping of 5bp
#' result <- map_bed(tbl_x, tbl_y, operation = list(score_mean = list(on = "score", func = mean)), 
#'                   min_overlap = 5, min_overlap_type = "bp")
#' head(result)
#' @references Manual page of `bedtools map`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/map.html}
#' @seealso [bedtorch::merge_bed()], [bedtorch::intersect_bed()] [GenomicRanges::findOverlaps()]
#' @export
map_bed <- function(x,
                    y,
                    operation,
                    max_gap = -1L,
                    min_overlap = 0L,
                    min_overlap_type = c("bp", "frac1", "frac2")) {
  UseMethod("map_bed")
}

#' @export
map_bed.bedtorch_table <- function(x,
                               y,
                               operation,
                               max_gap = -1L,
                               min_overlap = 0L,
                               min_overlap_type = c("bp", "frac1", "frac2")) {
  data <- x
  scaffold <- y
  stopifnot(is(data, "bedtorch_table"))
  stopifnot(is(scaffold, "bedtorch_table"))
  stopifnot(!is.null(operation))
  
  result <-
    intersect_bed.bedtorch_table(
      data,
      scaffold[, 1:3],
      mode = "wo",
      max_gap = max_gap,
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
  
  genome <- attr(result, "genome")
  op_names <- names(operation)
  result <- result[, {
    dt2 <- lapply(op_names, function(op_name) {
      func <- operation[[op_name]]$func
      col_name <- operation[[op_name]]$on
      func(.SD[[col_name]])
    })
    names(dt2) <- op_names
    dt2
  }, by = c("chrom", "start", "end")][scaffold[, 1:3], nomatch = 0]
  new_bedtorch_table(result, genome)
}


#' Subtract 
#' 
#' This function searches for intervals in `y` that overlaps with `x`, and
#' remove the overlapping portion of `x`.
#' @param x A `GRanges` object.
#' @param y A `GRanges` object.
#' @param min_overlap See [bedtorch::intersect_bed()].
#' @param min_overlap_type [bedtorch::intersect_bed()].
#' @return A `GRanges` after subtracting `y` from `x`
#' @examples
#' # Load BED tables
#' tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#' tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"), use_gr = FALSE)
#'
#' # Basic usage
#' result <- subtract_bed(tbl_x, tbl_y)
#' head(result)
#'
#' # Perform the mapping, requiring the minimum overlapping of 5bp
#' result <- subtract_bed(tbl_x, tbl_y, min_overlap = 5, min_overlap_type = "bp")
#' head(result)
#' @references Manual page of `bedtools subtract`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html}
#' @seealso [bedtorch::intersect_bed()]
#' @export
subtract_bed <- function(x, y, min_overlap = 1, min_overlap_type = c("bp", "frac1", "frac2")) {
  UseMethod("subtract_bed")
}

#' @export
subtract_bed.GRanges <- function(x, y, min_overlap = 1, min_overlap_type = c("bp", "frac1", "frac2")) {
  x <- as.bedtorch_table(x)
  y <- as.bedtorch_table(y)
  
  subtract_bed.bedtorch_table(x, y, min_overlap, min_overlap_type) %>%
    as.GenomicRanges()
}

#' @export
subtract_bed.bedtorch_table <- function(x, y, min_overlap = 1, min_overlap_type = c("bp", "frac1", "frac2")) {
  y <- merge_bed(y)
  overlap_scaffold <-
    intersect_bed(
      x[, 1:3],
      y[, 1:3],
      mode = "wo",
      min_overlap = min_overlap,
      min_overlap_type = min_overlap_type
    )
  result <- overlap_scaffold[x]
  
  helper <- function(.SD, chrom, start, end) {
    if (is.na(.SD[[1]][1])) {
      return(.SD)
    }
    
    start_s <- c(start, .SD[[3]])
    end_s <- pmax(c(.SD[[2]], end), start_s)
    
    valid_idx <- end_s > start_s
    start_s <- start_s[valid_idx]
    end_s <- end_s[valid_idx]
    cnt <- length(start_s)
    
    if (cnt == 0)
      return(NULL)
    
    rhs <- list(start_s = start_s, end_s = end_s)
    
    lhs <- lapply(seq_len(ncol(.SD) - 2), function(col_idx) {
      if (cnt > 1) {
        rep(.SD[[col_idx]][1], cnt)
      } else
        .SD[[col_idx]][1]
    })
    names(lhs) <- colnames(.SD)[seq_along(lhs)]
    c(lhs, rhs)
  }
  
  # Preprocess
  # Remove overlap
  result[, (7) := NULL]
  result[, `:=`(start_s = start, end_s = end)]
  
  result <- result[, {
    helper(.SD, chrom, start, end)
  }, by = c("chrom", "start", "end")]
  
  # Example
  # chrom start   end chrom.1 start.1 end.1 score1 score2 start_s end_s
  # 1:     1  1057  1146       1     842  1136     86    188    1136  1146
  # 2:     1  1130  1228       1     842  1136    116    185    1136  1151
  # 3:     1  1397  1486       1    1151  1413     86    193    1413  1418
  # 4:     1  1710  1799       1    1418  1796    110    199    1796  1799
  # 5:     1  1788  1895       1    1418  1796     96    176    1796  1805
  result[, 2:6 := NULL]
  setnames(result, c("start_s", "end_s"), c("start", "end"))
  result <-
    result %>% dplyr::relocate(c("start", "end"), .after = "chrom")
  setkey(result, "chrom", "start", "end")
  
  # Determine the genome attribute of the reuslt
  genome <- local({
    genome_list <-
      list(x, y) %>% map( ~ attr(., which = "genome")) %>% 
      discard(is.null) %>% unique()
    
    if (length(genome_list) == 2)
      # Dual genome, can't decide the genome of the result
      NULL
    else if (length(genome_list) == 0)
      NULL
    else
      genome_list[[1]]
  })
  new_bedtorch_table(result, genome)
}


#' Generate complement
#' 
#' This function returns all intervals in a genome that are not covered by at least one interval in the input.
#' @param x A `GRanges`.
#' @param genome Specify the reference genome for the BED file. `genome` can be
#'   a valid genome name in [GenomeInfoDb::Seqinfo()], e.g. `GRCh37`, or
#'   `hs37-1kg`, which is a genome shipped with this package, or any custom
#'   chromosome size files (local or remote). If `NULL`, will try to obtain such
#'   information from input data. Refer to [bedtorch::read_bed()].
#' @return A `GRanges` which is the complement of `x`.
#' @examples
#' # Load BED tables
#' tbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#'
#' # Basic usage
#' result <- complement_bed(tbl, "hs37-1kg")
#' head(result)
#' @references Manual page of `bedtools complement`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/complement.html}
#' @seealso [bedtorch::subtract_bed()]
#' @export
complement_bed <- function(x, genome = NULL) {
  UseMethod("complement_bed")
}

#' @export
complement_bed.GRanges <- function(x, genome = NULL) {
  x <- as.bedtorch_table(x)
  complement_bed.bedtorch_table(x, chrom_sizes) %>%
    as.GenomicRanges()
}

#' @export
complement_bed.bedtorch_table <- function(x, genome = NULL) {
  if (is.null(genome))
    genome <- attr(x, "genome")
  if (is.null(genome))
    stop("Must provide genome information")
  
  chrom_sizes <- load_chrom_sizes(genome)
  
  chrom_sizes <-
    chrom_sizes[, `:=`(start = 0L, end = as.integer(size))][, .(chrom, start, end)]
  setkey(chrom_sizes, "chrom", "start", "end")
  chrom_sizes <- new_bedtorch_table(chrom_sizes, genome = attr(x, "genome"))
  
  subtract_bed(chrom_sizes, x)
}


#' Shuffle intervals
#' 
#' This function randomly permutes genomic locations among a genome defined in a
#' genome file. Note: the permutation stays in the same chromosome.
#' Cross-chromosome shuffle is not supported at this point.
#' @param x A `GRanges`.
#' @param genome Specify the reference genome for the BED file. `genome` can be
#'   a valid genome name in [GenomeInfoDb::Seqinfo()], e.g. `GRCh37`, or
#'   `hs37-1kg`, which is a genome shipped with this package, or any custom
#'   chromosome size files (local or remote). If `NULL`, will try to obtain such
#'   information from input data. Refer to [bedtorch::read_bed()].
#' @param excluded_region A `GRanges` containing regions that you don't want
#'   to exclude from the shuffle. Shuffled intervals are guaranteed not to fall
#'   within such regions.
#' @param sort Logical value indicating whether to sort the shuffled features.
#'   When the input data set is large, sorting can be very expensive. If you
#'   don't need the output to be sorted, set `sort` to `FALSE`.
#' @param seed An integer seed for the random number generator. If `NULL`, the
#'   seed is randomly chosen.
#' @return A `GRanges` containing shuffled features.
#' @examples
#' # Load BED tables
#' tbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE, genome = "hs37-1kg")
#' excluded <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"), use_gr = FALSE)
#'
#' # Basic usage
#' result <- shuffle_bed(tbl)
#' head(result)
#' 
#' # Shuffle the data, but exclude certain regions. Plus, set the RNG seed to 1
#' result <- shuffle_bed(tbl, genome = "hs37-1kg", excluded_region = excluded, seed = 1)
#' head(result)
#' @references Manual page of `bedtools shuffle`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}
#' @export
shuffle_bed <-
  function(x,
           genome = NULL,
           excluded_region = NULL,
           sort = TRUE,
           seed = NULL) {
    UseMethod("shuffle_bed")
  }


#' @export
shuffle_bed.GRanges <-
  function(x,
           genome = NULL,
           excluded_region = NULL,
           sort = TRUE,
           seed = NULL) {
    x <- as.bedtorch_table(x)
    shuffle_bed.bedtorch_table(x, genome, excluded_region, sort, seed) %>%
      as.GenomicRanges()
  }


#' @export
shuffle_bed.bedtorch_table <-
  function(x,
           genome = NULL,
           excluded_region = NULL,
           sort = TRUE,
           seed = NULL) {
    if (length(unique(x$chrom)) > 1) {
      result <- rbindlist(unique(x$chrom) %>% map(function(chrom0) {
        shuffle_bed(x[chrom == chrom0], genome, excluded_region)
      }))
      new_bedtorch_table(result, genome = attr(x, "genome"))
    }

  # Single chrom
  chrom0 <- x$chrom[1]
  
  if (is.null(genome))
    genome <- attr(x, "genome")
  if (is.null(genome))
    stop("Must provide genome information")
  
  chrom_sizes <- load_chrom_sizes(genome)
  chrom_sizes <-
    chrom_sizes[, `:=`(start = 0L, end = as.integer(size))][chrom == chrom0, .(chrom, start, end)]
  setkey(chrom_sizes, "chrom", "start", "end")
  chrom_sizes <- new_bedtorch_table(chrom_sizes)
  chrom_sz <- chrom_sizes[, end - start]
  
  if (!is.null(excluded_region)) {
    excluded_region <- merge_bed(excluded_region[chrom == chrom0])
    available_width <- chrom_sz - excluded_region[, sum(end - start)]
    available_region <- subtract_bed(chrom_sizes, excluded_region)
    available_scaffold <- as.integer(c(0, available_region[, end - start] %>% cumsum()))
  } else {
    available_width <- chrom_sz
    available_region <- chrom_sizes
    available_scaffold <- c(0L, chrom_sz)
  }
  
  # Map the scaffold coordinate to the real genomic coordinate
  map2coord <- function(pos_scaffold, width) {
    # pos_scaffold falls into #block_id region
    pred1 <- available_scaffold > pos_scaffold
    pred2 <- available_scaffold >= pos_scaffold + width
    if (any(pred1 != pred2)) {
      # Start/end of the simulated fragment are in different blocks
      return(NA)
    }
    if (!any(pred1)) {
      # pos_scaffold is out of range
      return(NA)
    }
    
    block_id <- min(which(pred1)) - 1
    available_region[block_id, start + pos_scaffold - available_scaffold[block_id]]
  }
  
  
  # Generate the random intervals
  if (!is.null(seed))
    set.seed(seed = seed)

  x_scaffold <- x[, 1:3][, s1 := as.integer(NA)]
  while(TRUE) {
    # All fragments are shuffled
    if (nrow(x_scaffold[is.na(s1)]) == 0)
      break
    
    x_scaffold[is.na(s1), s1 := {
      random_start <- purrr::map2_int(sample.int(
        n = available_width,
        size = nrow(.SD),
        replace = TRUE
      ) - 1L,
      end - start,
      map2coord)
      random_start
    }][, e1 := s1 + (end - start)]
  }
  
  x_scaffold[, `:=`(start=NULL, end=NULL)]
  setnames(x_scaffold, c("s1", "e1"), c("start", "end"))
  result <- cbind(x_scaffold, x[, -1:-3])
  if (sort)
    setkey(result, "chrom", "start", "end")
  
  new_bedtorch_table(result, genome = attr(x_scaffold, "genome"))
}


#' Cluster intervals
#' 
#' This operation is similar to `bedtools cluster`. It reports each set of
#' overlapping or “book-ended” features in an interval file.
#' @param x A `GRanges`
#' @param max_dist Maximum distance between features allowed for features to be
#'   merged. Default is 0. That is, overlapping and/or book-ended features are
#'   merged.
#' @return A `GRanges`. Compared
#' @examples 
#' tbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#' 
#' # Basic usage
#' clustered <- cluster_bed(tbl)
#' head(clustered)
#' 
#' # Change the maximum distance allowed
#' clustered <- cluster_bed(tbl, max_dist = 10)
#' head(clustered)
#' @seealso [bedtorch::merge_bed()]
#' @references Manual page of `bedtools cluster`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/cluster.html}
#' @export
cluster_bed <- function(x, max_dist = 0) {
  UseMethod("cluster_bed")
}

#' @export
cluster_bed.GRanges <- function(x, max_dist = 0) {
  x <- as.bedtorch_table(x)
  cluster_bed.bedtorch_table(x, max_dist) %>%
    as.GenomicRanges()
}


#' @export
cluster_bed.bedtorch_table <- function(x, max_dist = 0) {
  stopifnot(max_dist >= 0)
  
  idx_colname <- .available_colname(x, "cluster")
  x_with_cluster <- x[, {
    cluster_id <-
      head(c(0, cumsum(
        shift(start, n = -1) > cummax(end) + max_dist
      )), n = -1) + 1L
    cluster_col <- list()
    cluster_col[[idx_colname]] <- cluster_id
    c(.SD, cluster_col)
  }, by = chrom]
  
  # Cluster id over chromosomes
  cluster_map <- unique(x_with_cluster[, c("chrom", ..idx_colname)])
  setnames(cluster_map, c("chrom", "cluster"))
  # This is the global cluster ID
  cluster_map[, g_cluster := seq_len(.N)]
  setkey(cluster_map, "chrom", "cluster")
  setkey(x_with_cluster, "chrom", "cluster")
  x_with_cluster = merge(x_with_cluster, cluster_map, all.x = TRUE)
  setkey(x_with_cluster, "chrom", "start", "end")
  x_with_cluster[, (idx_colname) := NULL]
  setnames(x_with_cluster, "g_cluster", "cluster")
  x_with_cluster[]
}


#' Make windows along genome
#' 
#' Make windows with the size of `window_size` along genome. Must provide
#' `genome` so that intervals are restricted by chromosome sizes.
#' @param window_size An integer value. Size of the windows in base pairs.
#' @param genome Provide the genome identifier, e.g. `GRCh37`, `hs37-1kg`, or
#'   path/URL to the chromosome size file. Refer to [GenomeInfoDb::Seqinfo()]
#'   and [bedtorch::read_bed()].
#' @param chrom A character vector. If provided, only generate windows for the
#'   specified chromosome(s).
#' @return A `GRanges` object containing the generated windows.
#' @examples 
#' # Generate 500kbp intervals for chr20 and chr22, GRCh38
#' result <- make_windows(window_size = 500e3L, genome = "GRCh38", chrom = c("20", "22"))
#' head(result)
#' 
#' # Generate 5Mbp intervals for all chromosomes of GRCh37 (default reference genome)
#' result <- make_windows(window_size = 5e6L, genome = "GRCh37")
#' head(result)
#' @export
make_windows <- function(window_size, genome, chrom = NULL) {
  stopifnot(!is.null(genome))
  
  # Get chrom_sizes from genome seqinfo
  genome <- get_seqinfo(genome)
  chrom_sizes <- genome %>% as.data.frame()
  chrom_sizes$chrom <- row.names(chrom_sizes)
  setDT(chrom_sizes)
  chrom_sizes <- chrom_sizes[, .(chrom, size = seqlengths)]
  
  chrom_list <- chrom
  if (!is.null(chrom_list))
    chrom_sizes <- chrom_sizes[chrom %in% chrom_list]
  
  stopifnot(nrow(chrom_sizes) > 0)
  
  
  result <- lapply(seq.int(nrow(chrom_sizes)), function(row_idx) {
    chrom <- chrom_sizes$chrom[row_idx]
    size <- chrom_sizes$size[row_idx]
    
    start_pos <-
      seq(1, size, by = window_size) %>% head(n = ceiling(size / window_size))
    end_pos <- pmin(start_pos + window_size - 1, size)
    GenomicRanges::GRanges(
      seqnames = chrom,
      ranges = IRanges::IRanges(start = start_pos, end = end_pos),
      seqinfo = genome
    )
  })
  do.call(c, result)
}


#' Randomly place intervals of fixed width across genome
#' 
#' This function is aware of chromosome size differences and intervals are
#' randomly generated in proportion to chromosome sizes. In another word, larger
#' chromosomes will have more randomly generated intervals.
#' @param n Number of intervals to generate.
#' @param interval_width An integer value of the width of the interval.
#' @param genome Specify the reference genome for the BED file. `genome` can be
#'   a valid genome name in [GenomeInfoDb::Seqinfo()], e.g. `GRCh37`, or
#'   `hs37-1kg`, which is a genome shipped with this package, or any custom
#'   chromosome size files (local or remote). If `NULL`, will try to obtain such
#'   information from input data. Refer to [bedtorch::read_bed()].
#' @param chrom A character vector. If provided, only generate random intervals
#'   for the specified chromosome(s).
#' @param seed An integer seed for the random number generator. If `NULL`, the
#'   seed is randomly chosen.
#' @param sort Logical value. Should the random BED data be sorted?
#' @return A `GRanges` of generated intervals.
#' @examples 
#' # Generate 100 500kbp-intervals for chr20 and chr22, GRCh38
#' result <- make_random_bed(100, 500e3L, genome = "GRCh38", chrom = c("20", "22"))
#' head(result)
#' 
#' # Generate 100 500bp unsorted intervals for all chromosomes of GRCh37 (default reference genome)
#' result <- make_random_bed(100, interval_width = 500, genome = "GRCh37", sort = FALSE)
#' head(result)
#' @references Manual page for `bedtools random`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/random.html}
#' @export
make_random_bed <- function(n, interval_width, genome, chrom = NULL, seed = NULL, sort = TRUE) {
  # BED for chrom_sizes
  if (is.null(genome))
    stop("Must provide genome information")
  
  chrom_sizes <- load_chrom_sizes(genome)
  
  chrom_list <- chrom
  if (!is.null(chrom_list))
    chrom_sizes <- chrom_sizes[chrom %in% chrom_list]
  
  # Generate n intervals and assign chrom
  sample_chrom <-
    chrom_sizes$chrom[sample.int(
      n = nrow(chrom_sizes),
      size = n,
      replace = TRUE,
      prob = chrom_sizes$size
    )]
  
  # Generate the random intervals
  if (!is.null(seed))
    set.seed(seed = seed)
  
  result <- chrom_sizes[, .(chrom = sample_chrom, v = 1)]
  result <- result[, {
    chrom0 <- chrom
    size <- chrom_sizes[chrom == chrom0]$size
    
    start = sample.int(n = size - interval_width + 1, size = nrow(.SD)) - 1L
    end = start + interval_width
    list(start = start, end = end)
  }, by = chrom]
  
  if (sort)
    setkey(result, "chrom", "start", "end")
  
  new_bedtorch_table(result, genome = genome)
}


#' Jaccard distance
#' 
#' Calculate Jaccard distance between two BED-like data sets
#' 
#' @param x A `GRanges`.
#' @param y A `GRanges`.
#' @param min_overlap The smallest overlapping region for two intervals to be
#'   considered as overlapping. Default is 1.
#' @param min_overlap_type A character value indicating how `min_overlap` is
#'   interpreted. `bp` means `min_overlap` is the number of base pairs. `frac1`
#'   means `min_overlap` is the minimum overlap required as a fraction of `x`.
#'   Similarly, `frac2` means `min_overlap` is the minimum overlap required as a
#'   fraction of `y`. Similar to `bedtools intersect`'s `-f` and `-F` arguments.
#' @return The Jaccard statistics.
#' @seealso [bedtorch::intersect_bed()] 
#' @examples 
#' # Load BED tables
#' x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"), use_gr = FALSE)
#' y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"), use_gr = FALSE)
#' 
#' # Basic usages
#' result <- jaccard_bed(x, y) 
#' head(result)
#' @references Manual page of `bedtools jaccard`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html}
#' @export
jaccard_bed <- function(x,
                        y,
                        min_overlap = 1,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  UseMethod("jaccard_bed")
}


#' @export
jaccard_bed.GRanges <- function(x,
                        y,
                        min_overlap = 1,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  x <- as.bedtorch_table(x)
  y <- as.bedtorch_table(y)
  
  jaccard_bed.bedtorch_table(x, y, min_overlap, min_overlap_type)
}


#' @export
jaccard_bed.bedtorch_table <- function(x,
                        y,
                        min_overlap = 1,
                        min_overlap_type = c("bp", "frac1", "frac2")) {
  x <- merge_bed(x)
  y <- merge_bed(y)
  
  genome <- local({
    genome_list <-
      list(x, y) %>% map( ~ attr(., which = "genome")) %>% 
      discard(is.null) %>% unique()
    
    if (length(genome_list) == 2)
      # Dual genome, can't decide the genome of the result
      NULL
    else if (length(genome_list) == 0)
      NULL
    else
      genome_list[[1]]
  })
  
  union_data <- rbindlist(list(x[, 1:3], y[, 1:3]))
  setkey(union_data, "chrom", "start", "end")
  union_data <- new_bedtorch_table(union_data, genome)
  union_data <- merge_bed(union_data)
  
  intersect_data <-
    intersect_bed(x, y, min_overlap = min_overlap, min_overlap_type = min_overlap_type)
  
  union_sum <- union_data[, sum(end - start)]
  intersect_sum <- intersect_data[, sum(end - start)]
  
  jaccard <- intersect_sum / union_sum
  
  data.table(
    intersection = intersect_sum,
    union = union_sum,
    jaccard = jaccard,
    n_intersections = nrow(intersect_data)
  )
}


#' Relative distance distribution
#' 
#' Calculate the relative distance distribution between two BED data tables.
#' @param x A `GRanges`.
#' @param y A `GRanges`.
#' @return The relative distance distribution.
#' @examples 
#' # Load BED tables
#' x <- read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"), use_gr = FALSE)
#' y <- read_bed(system.file("extdata", "example2_window.bed", package = "bedtorch"), use_gr = FALSE)
#' 
#' # Basic usages
#' reldist_bed(x, y) 
#' @references 
#' Manual page of `bedtools reldist`:
#'   \url{https://bedtools.readthedocs.io/en/latest/content/tools/reldist.html}
#'   
#' Exploring Massive, Genome Scale Datasets with the GenometriCorr Package.
#' Favorov A, Mularoni L, Cope LM, Medvedeva Y, Mironov AA, et al. (2012) PLoS
#' Comput Biol 8(5): e1002529.
#' \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529}{doi:10.1371/journal.pcbi.1002529}
#' @export
reldist_bed <- function(x, y) {
  UseMethod("reldist_bed")
}

#' @export
reldist_bed.GRanges <- function(x, y) {
  x <- as.bedtorch_table(x)
  y <- as.bedtorch_table(y)
  
  reldist_bed.bedtorch_table(x, y) %>%
    as.GenomicRanges()
}


#' @export
reldist_bed.bedtorch_table <- function(x, y) {
  helper <- function(x, y) {
    if (min(nrow(x), nrow(y)) == 0)
      return(NULL)

    # This helper function runs only for single chrom
    stopifnot(length(unique(x$chrom)) == 1)
    stopifnot(length(unique(y$chrom)) == 1)
    
    # Midpoints
    # x <- x[, 1:3][, midpoint := as.integer(round(start + (end - start - 1) / 2))][, type := "x"]
    # y <- y[, 1:3][, midpoint := as.integer(round(start + (end - start - 1) / 2))][, type := "y"]
    x <- x[, 1:3][, midpoint := as.integer(((start + end) / 2))][, type := "x"]
    y <- y[, 1:3][, midpoint := as.integer(((start + end) / 2))][, type := "y"]
    
    
    m <- rbindlist(list(x, y))
    setkey(m, "midpoint", "type")
    
    ptr1 <- NULL
    ptr2 <- NULL
    result <- rep(NA, nrow(m))
    
    seq_len(nrow(m)) %>% walk(function(idx) {
      if (m$type[idx] == "y") {
        if (!is.null(ptr1))
          ptr2 <<- c(ptr2, idx)
      } else {
        if (is.null(ptr1))
          ptr1 <<- idx
        else {
          if (is.null(ptr2))
            ptr1 <<- idx
          else {
            # Complete one cycle: have found the pattern x->y->x
            p1 <- m$midpoint[ptr1]
            p2 <- m$midpoint[ptr2]
            p3 <- m$midpoint[idx]
            
            # Update the result vector
            if (any(is.na(c(p1, p2, p3))))
              browser()
            if (any(is.null(c(p1, p2, p3))))
              browser()
            reldist <- pmin(p2 - p1, p3 - p2) / (p3 - p1)
            if (any(is.infinite(reldist)))
              browser()
            result[ptr2] <<- reldist
            
            ptr1 <<- idx
            ptr2 <<- NULL
          }
        }
      }
    })
    
    m <- m[, reldist := result][!is.na(result)]
    m$reldist
  }
  
  result <- unique(c(as.character(x$chrom), as.character(y$chrom))) %>%
    map(function(v) {
      helper(x[chrom == v], y[chrom ==v])
    }) %>%
    unlist()

  h <- hist(result, breaks = seq(0, 0.5, by = 0.01))
  data.table(
    reldist = head(h$breaks, n = -1),
    count = h$counts,
    total = length(result)
  )[, freq := count / total][count != 0]
}