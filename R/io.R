# Check if a file name is a gzip file
# 
# is_gzip("foo.gz")
# is_gzip("http://foo.bar/example.bed.gz")
# is_gzip(c("foo.gz", "bar.gz"))
is_gzip <- function(file_path) {
  str_split(file_path,
            pattern = "\\.") %>% map_chr( ~ tail(., n = 1)) %in% c("gz", "bgz")
}

is_remote <- function(file_path) {
  str_detect(file_path, pattern = "^[^:]+://")
}


# Make sure all ranges are valid
normalize_tabix_range <- function(range) {
  stopifnot(all(seqminer::isTabixRange(range)))
  
  # If some of the ranges are in the form chr1 (while chromosome), change it to
  # chr1:1-2000000000, because the underlying tabix library does not allow
  # single-chromosome range notation
  range_splits <- str_split(range, pattern = ":")
  single_chrom <- range_splits %>% map_lgl(~ length(.) == 1)
  
  # R's integer type is 32-bit, so the largest range is approx. 2 billion
  range[single_chrom] <-
    range_splits[single_chrom] %>%
    map_chr(function(v) {
      str_interp("${v[1]}:1-2100000000")
    })
  range
}


# Set `dt`'s column names, types, and set the index
# Modifies `dt` in-place.
post_process_table <- function(dt) {
  default_colnames <- all(str_detect(colnames(dt), pattern = "V[0-9]+"))
  
  # First 3 columns
  data.table::setnames(dt, old = 1:3, new = c("chrom", "start", "end"))
  
  # Fourth columns: score if it is numeric, feature otherwise
  if (default_colnames && ncol(dt) >= 4) {
    if (is.numeric(dt[[4]]))
      data.table::setnames(dt, 4, "score")
    else
      data.table::setnames(dt, 4, "feature")
  }
  
  v_chrom <- as.character(dt$chrom)
  chrom_list <- str_sort(unique(v_chrom), numeric = TRUE)
  dt[, chrom := factor(v_chrom, levels = chrom_list)]
  dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  data.table::setkey(dt, "chrom", "start", "end")
  dt[]
}


#' Check if required binaries exist in PATH
#'
#' Some of cragr's functionalities require certain third-party binaries to exist
#' in PATH. For example, selectively loading cfDNA fragments from a certain
#' region under the hood invokes tabix to complete the task. This function checks
#' whether certain binaries exist.
#' @param binaries Character vector indicating what binaries to check. Default
#'   is tabix.
#' @param verbose Logical value indicating whether to output detailed messages.
#'   Default is FALSE.
#' @param stop_on_fail Stop the script completely if the result is `FALSE`.
#' @return Logical value. TRUE if all required binaries can be found in PATH.
#' @export
check_binaries <- function(binaries = "tabix", verbose = FALSE, stop_on_fail = FALSE) {
  check_results <- binaries %>% purrr::map_lgl(function(binary) {
    # check if binary is in path
    if ("" == Sys.which(binary)) {
      if (verbose) {
        message(str_interp("Checking ${binary}: FAIL"))
      }
      return(FALSE)
    }
    else {
      if (verbose) {
        message(str_interp("Checking ${binary}: PASS"))
      }
      return(TRUE)
    }
  })
  
  check_results <- all(check_results)
  if (stop_on_fail)
    assertthat::assert_that(check_results, msg = paste0("Failed in locating ", paste(binaries, collapse = ", ")))
  else
    check_results
}


filter_by_region <- function(dt, range) {
  stopifnot(length(range) == 1)
  
  if (str_detect(range, pattern = "^[^:]+$"))
    return(dt[chrom == range])
  
  m <- str_match(range, pattern = "([^:]+):([0-9]+)-([0-9]+)")
  if (is.na(m[1, 1]))
    stop(str_interp("Invalid range: ${range}"))
  else {
    c1 <- m[1, 2]
    s1 <- as.integer(m[1, 3]) - 1L
    e1 <- as.integer(m[1, 4])
    return(dt[chrom == c1 & ((end > s1 & end <= e1) | (start < e1 & start >= s1))])
  }
}


# Load a gzipped and indexed BED-like file
read_tabix_bed <- function(file_path, range, use_htslib = FALSE) {
  # Make sure the genomic ranges are valid
  range <- normalize_tabix_range(range)

  na_strings <- c("NA", "na", "NaN", "nan", ".", "")
  filter_flag <- FALSE
  if (!use_htslib) {
    # Check if we can use tabix to extract the reads
    tabix_exist <- check_binaries(binaries = "tabix", verbose = FALSE, stop_on_fail = TRUE)
    index_exist <- file.exists(paste0(file_path, ".tbi"))
    use_tabix <- tabix_exist && index_exist
    
    if (use_tabix) {
      tabix_regions <- paste0(as.character(range), collapse = " ")
      tabix_cmd <- str_interp("tabix ${file_path} ${tabix_regions}")
      dt <-
        data.table::fread(cmd = tabix_cmd, na.strings = na_strings)
    } else {
      if (endsWith(file_path, ".gz")) {
        warning("Either tabix or index does not exist. Resource to sequential scanning, which may negatively impact the performance.")
      }
      dt <- data.table::fread(file = file_path, na.strings = na_strings)
      # Postpone the filtering to after assigning column names
      filter_flag <- TRUE
    }
  } else {
    dt <- seqminer::tabix.read.table(file_path, tabixRange = range)
    data.table::setDT(dt)
  }
  
  if (nrow(dt) == 0)
    return(dt)
  
  post_process_table(dt)
  
  if (filter_flag) {
    dt <- filter_by_region(dt, range)
  }
  
  dt
}


#' Load a BED-format file
#'
#' This function loads the input file as a `data.table` object. The file can be
#' either local or remote, and can be either plain text or gzip-compressed.
#' Furthermore, this function supports range-loading by providing a genomic
#' range in the following syntax: "chr1:1-100".
#' 
#' Note: for loading remote data files, currently this function depends on
#' tabix.c 0.2.5, which doesn't not support HTTPS protocol. In the next step, I
#' plan to turn to htslib, and the this function can load remote data files
#' through HTTPS.
#' 
#' @param file_path Path to the data file. It can be either a local file, or a remote URL.
#' @param range A genomic range character vector. Must follow standard genomic
#'   range notation format, e.g. chr1:1001-2000
#' @param ... Other arguments to be passed to [data.table::fread()].
#' @seealso [data.table::fread()]
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' head(bedtbl)
#' 
#' bedtbl <- read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
#'                    range = "1:3001-4000")
#' head(bedtbl)
#' 
#' # Does not work. Currently HTTPS is not supported
#' # bedtbl <- read_bed("https://yourdomain.com/example-02.bedGraph.gz", 
#' #                    range = "1:1000-2000")
#' @export
read_bed <- function(file_path, range = NULL, ...) {
  if (is.null(range)) {
    # Load directly
    na_strings <- c("NA", "na", "NaN", "nan", ".", "")
    dt <- fread(file_path, sep = "\t", na.strings = na_strings, ...)
    post_process_table(dt)
  } else {
    stopifnot(length(range) == 1)
    
    if (is_gzip(file_path)) {
      read_tabix_bed(file_path, range)
    } else {
      if (is_remote(file_path))
        stop("range filtering is not available for remote BGZIP files")
      else {
        # Load directly
        na_strings <- c("NA", "na", "NaN", "nan", ".", "")
        dt <- fread(file_path, sep = "\t", na.strings = na_strings, ...)
        post_process_table(dt)
        
        filter_by_region(dt, range)
      }
    }
  }
}


#' Write a `data.table` to file
#' 
#' Can be a plain text file or a BGZIP file.
#' 
#' @param x A `data.table` object to write.
#' @param file_path Path of the output file. If the extension name is `.gz` or
#'   `.bgz`, the output will be compressed in BGZIP format.
#' @param tabix_index If `TRUE`, and `file_path` indicates a gzip file, will
#'   create create the tabix index file.
#' @param ... Other arguments passed to methods. Compliant with `data.table::fwrite`.
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' 
#' # Write data to uncompressed file
#' write_bed(bedtbl,  tempfile(fileext = ".bed"))
#' 
#' # Write data to file and create tabix index
#' write_bed(bedtbl, tempfile(fileext = ".bed.gz"), tabix_index = TRUE)
#' @export
write_bed <- function(x, file_path, tabix_index = TRUE, ...) {
  compressed <- is_gzip(file_path)
  
  setnames(x, 1, "#chrom")
  on.exit(setnames(x, "#chrom", "chrom"), add = TRUE)
  
  if (compressed) {
    # Since we need to write the data table to disk as a temporary file, it's
    # important to operate by batches, i.e. in each batch, process rows no more
    # than batch_size
    batch_size <- 1e6L
    batch_plan <- seq(from = 1, to = nrow(x), by = batch_size)
    if (tail(batch_plan, n = 1) != nrow(x))
      batch_plan <- c(batch_plan, nrow(x))
    
    1:(length(batch_plan) - 1) %>%
      walk(function(batch_idx) {
        temp_txt <- tempfile(fileext = ".tsv")
        temp_gz <- tempfile(fileext = ".gz")
        on.exit(unlink(c(temp_txt, temp_gz)), add = TRUE)
        
        if (batch_idx < length(batch_plan) - 1) {
          # Not the last batch
          batch_data <- x[batch_plan[batch_idx]:(batch_plan[batch_idx + 1] - 1)]
        } else {
          batch_data <- x[batch_plan[batch_idx]:batch_plan[batch_idx + 1]]
        }
        
        data.table::fwrite(batch_data,
                           file = temp_txt,
                           quote = FALSE,
                           sep = "\t",
                           na = ".",
                           # Output column names only for batch #1
                           col.names = (batch_idx == 1),
                           ...)
        Rsamtools::bgzip(temp_txt, dest = temp_gz)
        unlink(temp_txt)
        
        if (batch_idx == 1)
          system(str_interp("cat ${temp_gz} > ${file_path}"))
        else
          system(str_interp("cat ${temp_gz} >> ${file_path}"))
      })
    
    if (tabix_index) {
      Rsamtools::indexTabix(file_path, format = "bed", zeroBased = TRUE)
    }
    ret <- TRUE
  } else {
    data.table::fwrite(x,
                       file = file_path,
                       quote = FALSE,
                       sep = "\t",
                       na = ".",
                       ...)
  }
}
