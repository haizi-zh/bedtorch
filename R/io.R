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


isTabixRange <- function(range) {
  isValid <- function(x) {
    ret = strsplit(x = x, split = ":")[[1]]
    if (length(ret) == 1) { # e.g. "chr1"
      return(TRUE)
    }
    if (length(ret) != 2) {
      return(FALSE)
    }
    chrom = ret[1]
    if (nchar(chrom) == 0) {
      return (FALSE)
    }
    ret = strsplit(x = ret[2], split = "-")[[1]]
    if (length(ret) == 2) {
      beg = suppressWarnings(as.integer(ret[1]))
      end = suppressWarnings(as.integer(ret[2]))
      if (is.na(beg) || is.na(end) || beg > end) {
        return(FALSE)
      }
    } else if (length(ret) == 1) {
      beg = suppressWarnings(as.integer(ret[1]))
      if (is.na(beg)) {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
    return(TRUE)
  }
  ranges <- unlist(strsplit(x = range, split = ","))
  sapply(ranges, isValid)
}


# Make sure all ranges are valid
normalize_tabix_range <- function(range) {
  stopifnot(all(isTabixRange(range)))
  
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


parse_range <- function(range) {
  # Make sure the genomic ranges are valid
  range <- normalize_tabix_range(range)
  
  range <- str_match(range, 
                     pattern = "^([^:]+):([0-9]+)-([0-9]+)")
  range_bed <- as.data.table(range)[, 2:4]
  setnames(range_bed, c("chrom", "start", "end"))
  range_bed[, `:=`(start = as.integer(start) - 1L, end = as.integer(end))]
  post_process_table(range_bed)
  new_bedtorch_table(range_bed) %>% merge_bed()
}


# Load a gzipped and indexed BED-like file
read_tabix_bed <- function(file_path, range, index_path = NULL, download_index = FALSE, ...) {
  # Get UCSC-style range strings
  range %<>%
    parse_range %>%
    pmap_chr(function(chrom, start, end)
      str_interp("${chrom}:${start + 1L}-${end}"))
  
  # Check index existence
  if (!is.null(index_path)) {
    if (is_remote(index_path))
      index_exists <- RCurl::url.exists(index_path)
    else
      index_exists <- file.exists(index_path)
    
    if (!index_exists) {
      warning(str_interp("Cannot find tabix index ${index_path}"))
      index_path <- NULL
    }
  } else if (is_remote(file_path)) {
    index_path <- paste0(file_path, ".tbi")
    index_exists <- RCurl::url.exists(index_path)
    if (!index_exists)
      stop(str_interp("Cannot find remote tabix index ${index_path}"))
  }
 
  tempbed <- tempfile(fileext = ".bed")
  on.exit(unlink(tempbed), add = TRUE)
  c_read_tabix_table(
    file_path,
    range,
    output_file = tempbed,
    index_path = if (is.null(index_path)) "" else index_path,
    download_index = download_index
  )
  # user_gr should be FALSE, since here we need a data.table object, which later
  # will be converted to GenomicRanges, if instructed
  
  # Load directly
  na_strings <- c("NA", "na", "NaN", "nan", ".", "")
  dt <-
    fread(tempbed, sep = "\t", na.strings = na_strings, ...)
  dt <- post_process_table(dt)
  dt
}


# Get a `GenomeInfoDb::Seqinfo` object
# 
# @param genome A canonical genome name, e.g. GRCh37, or a custom genome name
#   recognized by bedtorch, e.g. hs37-1kg, or a path/URL to a chromosome size
#   file. If `NULL`, return `NULL`.
# @param genome_name Optional character vector. Only works when `genome` is a
#   file path/URL, and specify the name of the genome. If `NULL`, the name is
#   guessed from the file name.
get_seqinfo <- function(genome, genome_name = NULL) {
  if (is.null(genome))
    return(NULL)
  
  if (file.exists(genome) || RCurl::url.exists(genome)) {
    # genome is a chromosome size file (local or remote)
    chrom_sizes <-
      fread(
        genome,
        col.names = c("chrom", "size")
      )
    
    if (is.null(genome_name)) {
      guess_genome_name <- function(genome) {
        # Guess genome name from file name
        genome_name <- basename(genome)
        # ends with chrom.sizes, sizes, or chromosome.sizes
        pattern <- "^(.+?)(\\.chrom(osome)?)?(\\.sizes)$"
        m <- str_match(genome_name, pattern = pattern)[1, 2]
        if (!is.na(m))
          m
        else
          genome_name
      }
      genome_name <- guess_genome_name(genome)
    }
    
    GenomeInfoDb::Seqinfo(
      seqnames = chrom_sizes$chrom,
      seqlengths = chrom_sizes$size,
      isCircular = rep(FALSE, nrow(chrom_sizes)),
      genome = genome_name %||% basename(genome)
    )
  } else if (genome == "hs37-1kg") {
    chrom_sizes <-
      fread(
        system.file("extdata", "human_g1k_v37.chrom.sizes", package = "bedtorch"),
        col.names = c("chrom", "size")
      )
    GenomeInfoDb::Seqinfo(
      seqnames = chrom_sizes$chrom,
      seqlengths = chrom_sizes$size,
      isCircular = rep(FALSE, nrow(chrom_sizes)),
      genome = genome
    )
  } else
    GenomeInfoDb::Seqinfo(genome = genome)
}


# Detect file format from a URL
detect_remote_file_type <- function(url) {
  # Download the first 4kb
  url_conn <- curl::curl(url)
  open(url_conn, "rb", blocking = FALSE)
  on.exit(close(url_conn), add = TRUE)
  
  total_bytes <- 0
  block_size <- as.integer(4 * 1024)
  total_buf <- list()
  while (isIncomplete(url_conn)) {
    buf <- readBin(url_conn, raw(), block_size)
    total_bytes <- total_bytes + length(buf)
    
    if (length(buf) == 0)
      next
    
    total_buf <- c(total_buf, list(buf))
    if (total_bytes >= block_size)
      break
  }
  total_buf <- unlist(total_buf)
  
  # Write to a temporary file and check the file type
  file_header <- tempfile()
  on.exit(unlink(file_header), add = TRUE)
  
  writeBin(total_buf, con = file_header)
  conn <- file(file_header)
  file_type <- summary(conn)$class
  close(conn)
  
  file_type
}


# Detect format of a local file
detect_local_file_type <- function(file_path) {
  conn <- file(file_path)
  file_type <- summary(conn)$class
  close(conn)
  
  file_type
}


# Read a URL without range seeking
read_bed_remote_full <- function(url, ...) {
  # Some URLs don't look like a file name. In this case, it's hard to detect the
  # file type solely from the URL Alternatively, we can download it to a
  # temporary file and start from there
  
  # Try to extract the extension name
  ext <- str_match(url, pattern = "\\.([^/\\.]+)$")[1, 2]
  if (is.na(ext))
    ext <- ""
  
  temp_downloaded <- tempfile(fileext = ext)
  on.exit(unlink(temp_downloaded), add = TRUE)
  
  curl::curl_download(
    url,
    temp_downloaded,
    mode = "wb",
    quiet = !getOption("datatable.showProgress", interactive())
  )
  
  na_strings <- c("NA", "na", "NaN", "nan", ".", "")
  
  # Guess file type
  conn <- file(temp_downloaded)
  file_type <- summary(conn)$class
  close(conn)
  if (file_type == "gzfile") {
    if (!endsWith(temp_downloaded, suffix = ".gz")) {
      # Make sure the file extension name is correct
      temp_download_new <- paste0(temp_downloaded, ".gz")
      file.rename(temp_downloaded, temp_download_new)
      on.exit(unlink(temp_download_new), add = TRUE)
      temp_downloaded <- temp_download_new
    }
    fread(file = temp_downloaded,
          na.strings = na_strings,
          sep = "\t",
          ...)
  } else if (file_type == "bzfile") {
    if (!endsWith(temp_downloaded, suffix = ".bz2")) {
      # Make sure the file extension name is correct
      temp_download_new <- paste0(temp_downloaded, ".bz2")
      file.rename(temp_downloaded, temp_download_new)
      on.exit(unlink(temp_download_new), add = TRUE)
      temp_downloaded <- temp_download_new
    }
    fread(file = temp_downloaded,
          na.strings = na_strings,
          sep = "\t",
          ...)
  } else if (file_type == "xzfile") {
    # LZMA or XZ
    fread(
      cmd = paste0("xzcat < ", temp_downloaded),
      na.strings = na_strings,
      sep = "\t",
      ...
    )
  } else {
    fread(file = temp_downloaded,
          na.strings = na_strings,
          sep = "\t",
          ...)
  }
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
#' @param compression Indicate the compression type. If `detect`, this function
#'   will try to guess from `file_path`.
#' @param tabix_index A character value indicating the location of the tabix
#'   index file. Can be either local or remote. If `NULL`, it will be derived
#'   from `file_path`.
#' @param download_index Whether to download (cache) the tabix index at current
#'   directory.
#' @param genome Specify the reference genome for the BED file. `genome` can be
#'   a valid genome name in [GenomeInfoDb::Seqinfo()], e.g. `GRCh37`, or
#'   `hs37-1kg`, which is a genome shipped with this package, or any custom
#'   chromosome size files (local or remote). Here is a good resource for such
#'   files: \url{https://github.com/igvteam/igv/tree/master/genomes/sizes}.
#' @param use_gr If `TRUE`, will read the data as a `GenomicRanges` object,
#'   otherwise a `data.table` object. Generally, we recommend using
#'   `GenomicRanges`.
#' @param ... Other arguments to be passed to [data.table::fread()].
#' @seealso [data.table::fread()]
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' head(bedtbl)
#' 
#' # Basic usage
#' bedtbl <- read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
#'                   range = "1:3001-4000")
#' head(bedtbl)
#' 
#' # Specify the reference genome
#' head(read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
#'               range = "1:3001-4000",
#'               genome = "hs37-1kg"))
#' 
#' head(read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
#'               range = "1:3001-4000",
#'               genome = "GRCh37"))
#' 
#' head(read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
#'               range = "1:3001-4000",
#'               genome = "https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/1kg_v37.chrom.sizes"))
#' 
#' # Load remote BGZIP files with tabix index specified
#' head(read_bed("https://git.io/JYATB", range = "22:20000001-30000001", tabix_index = "https://git.io/JYAkT"))
#' @export
read_bed <-
  function(file_path,
           range = NULL,
           compression = c("detect", "bgzip", "text", "other"),
           tabix_index = NULL,
           download_index = FALSE,
           genome = NULL,
           use_gr = TRUE,
           ...) {
    compression <- match.arg(compression)
    if (compression == "detect") {
      file_type <- if (is_remote(file_path)) 
        detect_remote_file_type(file_path)
      else
        detect_local_file_type(file_path)
      
      if (file_type == "gzfile")
        compression <- "bgzip"
      else
        compression <- "other"
    }
    
    if (is.null(range)) {
      # Load directly
      na_strings <- c("NA", "na", "NaN", "nan", ".", "")
      
      if (is_remote(file_path))
        dt <- read_bed_remote_full(file_path)
      else
        dt <- fread(file_path, sep = "\t", na.strings = na_strings, ...)
      dt <- post_process_table(dt)
    } else {
      stopifnot(length(range) == 1)
      
      if (compression == "bgzip") {
        dt <- read_tabix_bed(file_path,
                             range,
                             index_path = tabix_index,
                             download_index = download_index, ...)
      } else {
        if (is_remote(file_path))
          stop("Remote range filtering is only available for BGZIP files")
        else {
          # Load directly
          na_strings <- c("NA", "na", "NaN", "nan", ".", "")
          dt <-
            fread(file_path, sep = "\t", na.strings = na_strings, ...)
          post_process_table(dt)
          
          dt <- filter_by_region(dt, range)
        }
      }
    }
    
    dt <- new_bedtorch_table(dt, genome = genome)
    
    if (use_gr)
      as.GenomicRanges(dt)
    else
      dt
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
#' @param batch_size An positive integer. Write file in batches. During each
#'   batch, only write up to `batch_size` lines.
#' @param comments A character vector, which will be written to the top of the
#'   file as header lines.
#' @param ... Other arguments passed to methods. Compliant with `data.table::fwrite`.
#' @examples 
#' bedtbl <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
#' 
#' # Write data to uncompressed file
#' write_bed(bedtbl,  tempfile(fileext = ".bed"))
#' 
#' # Write data to file and create tabix index
#' write_bed(bedtbl, tempfile(fileext = ".bed.gz"), tabix_index = TRUE)
#' 
#' # Write data to uncompressed file, with header lines
#' write_bed(bedtbl,  tempfile(fileext = ".bed"), comments = c("Author: X", "Date: N/A"))
#' @export
write_bed <- function(x, file_path, tabix_index = TRUE, batch_size = NULL, comments = NULL, ...) {
  UseMethod("write_bed")
}


#' @export
write_bed.data.table <- function(x, file_path, tabix_index = TRUE, batch_size = NULL, comments = NULL, ...) {
  write_bed_core(x, file_path, tabix_index, batch_size, comments = comments, ...)
}


#' @export
write_bed.GRanges <- function(x, file_path, tabix_index = TRUE, batch_size = NULL, comments = NULL, ...) {
  dt <- data.table::as.data.table(x)
  dt[, `:=`(start = as.integer(start - 1L), width = NULL, strand = NULL)]
  write_bed_core(dt, file_path, tabix_index, batch_size, comments = comments, ...)
}


write_bed_core <- function(x, file_path, tabix_index = TRUE, batch_size = NULL, comments = NULL, ...) {
  compressed <- is_gzip(file_path)
  
  if (is(x, "GRanges")) {
    x <- data.table::as.data.table(x)
    x[, start := as.integer(start - 1L)]
  }
  
  setnames(x, 1, "#chrom")
  on.exit(setnames(x, "#chrom", "chrom"), add = TRUE)
  
  if (compressed) {
    # Since we need to write the data table to disk as a temporary file, it's
    # important to operate by batches, i.e. in each batch, process rows no more
    # than batch_size
    if (is.null(batch_size))
      batch_size <- nrow(x)
    batch_plan <- seq(from = 1, to = nrow(x), by = batch_size)
    if (tail(batch_plan, n = 1) != nrow(x))
      batch_plan <- c(batch_plan, nrow(x))
    
    1:(length(batch_plan) - 1) %>%
      walk(function(batch_idx) {
        temp_txt <- tempfile(fileext = ".tsv")
        # temp_gz <- tempfile(fileext = ".gz")
        on.exit(unlink(temp_txt), add = TRUE)
        
        if (batch_idx < length(batch_plan) - 1) {
          # Not the last batch
          batch_data <- x[batch_plan[batch_idx]:(batch_plan[batch_idx + 1] - 1)]
        } else {
          batch_data <- x[batch_plan[batch_idx]:batch_plan[batch_idx + 1]]
        }
        
        # Process comment lines
        if (!is.null(comments) && batch_idx == 1) {
          conn <- file(temp_txt)
          comments %>% 
            map_chr(function(x) paste0("#", x)) %>%
            writeLines(temp_txt)
          close(conn)
        }
        
        data.table::fwrite(batch_data,
                           file = temp_txt,
                           append = !is.null(comments) && batch_idx == 1,
                           quote = FALSE,
                           sep = "\t",
                           na = ".",
                           # Output column names only for batch #1
                           col.names = (batch_idx == 1),
                           ...)
        bgzip(temp_txt, output_file_path = file_path, append = (batch_idx > 1))
      })
    
    if (tabix_index) {
      build_tabix_index(file_path)
    }
  } else {
    # Process comment lines
    if (!is.null(comments)) {
      conn <- file(file_path)
      comments %>% 
        map_chr(function(x) paste0("#", x)) %>%
        writeLines(conn)
      close(conn)
    }
    
    data.table::fwrite(x,
                       file = file_path,
                       col.names = TRUE,
                       append = !is.null(comments),
                       quote = FALSE,
                       sep = "\t",
                       na = ".",
                       ...)
  }
}
