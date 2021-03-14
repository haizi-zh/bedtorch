# Home-made essentials tools for bioinformatics

#' @import readr
#' @import stringr
#' @import dplyr
#' @importFrom purrr modify_at map_chr
NULL


#' Check whether a url is remote or local
is_remote_url <- function(url) {
  schema <- str_match(url, "([a-z0-9]+)://.+")[, 2]
  
  sapply(schema, function(s) {
    if (is.na(s))
      FALSE
    else {
      if (!is.na(match(s, c("ftp", "sftp", "http", "https"))))
        TRUE
      else if (s == "file")
        FALSE
      else
        stop(paste0("Unknown schema: ", s))
    }
  })
}


#' Check if a file name is a gzip file
#' 
is_gzip_file <- function(file_path) {
  str_split(file_path,
            pattern = "\\.") %>% map_chr( ~ tail(., n = 1)) == "gz"
}


#' Load GENBED using tabix
#' 
#' @param file_path Path to the data file. Should be compressed by BGZIP,
#' locally present, and come along with a paired tabix index file
load_genbed_tabix <- function(file_path,
                              region,
                              col_names = NULL,
                              na = c("", "NA"),
                              col_types = NULL,
                              is_bedgraph = is_bedgraph,
                              ...) {
  genbed <-
    bedr::tabix(
      region = region,
      file.name = file_path,
      check.chr = FALSE,
      verbose = FALSE
    ) %>% as_tibble()
  
  if (!is.null(col_names)) {
    colnames(genbed) <- col_names
  } else if (all(colnames(genbed) == paste0("V", seq_along(colnames(genbed))))) {
    # Columns are default values: V1, V2, ...
    colnames(genbed) <- paste0("X", seq_along(colnames(genbed)))
    colnames(genbed)[1:3] <- c("chrom", "start", "end")
  }
  
  if (is_bedgraph) {
    # Try to convert the 4th column to numeric
    # In the original data table, identify rows with NA scores. Then identify
    # this kind of rows in converted score. If the indices match, the string-to-number
    # conversion is successful.
    na_idx <- genbed[[4]] %in% na
    score_dbl <-  as.numeric(genbed[[4]][!na_idx])
    if (all(!is.na(score_dbl))) {
      genbed[[4]] <- NA
      genbed[[4]][!na_idx] <- score_dbl
    }
  }
  
  genbed %>%
    modify_at(1, as.character) %>%
    modify_at(2:3, as.integer)
}


#' Load a data file in GENBED format
#'
#' This function loads a genbed file, which is simply a BED-like tsv file,
#' with a header line (optional) containing column names, and the types of
#' the first three columns are guaranteed to be "cii"
#' 
#' @param file_path Path or URL to the data file. Can be a local path, file://, 
#' http://, https://, ftp:// and sftp://
#' @param col_names can be NULL (infer from file), or a character vector
#' (user-provided column names). In the infer-from-file scenario, if there is 
#' a column header, we will have the column names. However, if the column header
#' is missing, then we will have the following column names: chrom, start, end,
#' X4, X5, ...
#' @param na Character vector of strings to interpret as missing values.
#' @param col_types A character string represents column names. Refer to
#' \code{col_type} in readr. If NULL, then the types for the first three columns
#' are character, integer, integer, and the rest are guessed from the data.
#' @param is_bedgraph A logic value indicating whether the input is BED or
#' BEDGRAPH. The difference between these two is that, the 4th column is character
#' for BED, and numeric for BEDGRAPH
#' @export
load_genbed <-
  function(file_path,
           region = NULL,
           col_names = NULL,
           na = c("", "NA"),
           col_types = NULL,
           is_bedgraph = FALSE,
           ...) {
    stopifnot(length(file_path) == 1)
    
    get_col_names <- function(file_path) {
      con <- file(file_path, "r")
      last_comment_line <- NULL
      first_data_line <- NULL
      n_fields <- 0
      while (TRUE) {
        line = str_trim(readLines(con, n = 1))
        
        if (length(line) == 0)
          break
        
        if (startsWith(line, "#")) {
          last_comment_line <- line
        }
        else {
          first_data_line <- line
          break
        }
      }
      close(con)
      
      # Number of fields
      if (is.null(first_data_line)) {
        # Empty data frame
        return(NULL)
      } else {
        n_fields <-
          length(str_split(string = first_data_line, pattern = "\t")[[1]])
      }
      
      stopifnot(n_fields >= 3)
      
      # Try extracting column names
      if (is.null(last_comment_line)) {
        col_names <- paste("X", 1:n_fields, sep = "")
        col_names[1:3] <- c("chrom", "start", "end")
      } else {
        col_names <- str_replace(string = last_comment_line,
                                 pattern = "^[# \t]+",
                                 replacement = "") %>%
          str_split(pattern = "\t") %>% .[[1]]
      }
      
      return(col_names)
    }
    
    tryCatch({
      temp_dir <- tempdir(check = TRUE)
      temp_path <- NULL
      temp_index_path <- NULL
      is_remote <- is_remote_url(file_path)
      
      # For remote files (either the data or the index), download them to 
      # the temporary dir 
      if (is_remote) {
        url <- file_path
        temp_path <-
          paste0(temp_dir, "/", tail(str_split(url, "/")[[1]], n = 1))
        download.file(url = url, destfile = temp_path)
        file_path <- temp_path
        
        # Attempt to download the tabix tbi file
        if (endsWith(url, ".gz") && !is.null(region)) {
          index_url <- paste0(url, ".tbi")
          if (RCurl::url.exists(index_url)) {
            temp_index_path <- paste0(temp_path, ".tbi")
            download.file(url = index_url, destfile = temp_index_path)
          }
        }
      }
      
      # Try identify column names from the header
      if (is.null(col_names)) {
        col_names <- get_col_names(file_path = file_path)
      }
      
      if (is.null(col_types)) {
        col_types <-
          paste(c("cii", rep_len("?", length(col_names) - 3)), collapse = "")
      }
      
      if (is.null(region)) {
        genbed <- read_tsv(
          file = file_path,
          comment = "#",
          col_names = col_names,
          col_types = col_types,
          na = na,
          ...
        )
      } else{
        genbed <-
          load_genbed_tabix(
            file_path = file_path,
            region = region,
            na = na,
            col_names = col_names,
            is_bedgraph = is_bedgraph,
            ...
          )
      }
      
      return(genbed)
    }, finally = {
      unlink(c(temp_path, temp_index_path, temp_dir))
    })
  }


# Dump the data frame as GENBED files
# Unlike
#
#' @export
write_genbed <-
  function(df,
           file_path,
           header = NULL,
           append = FALSE,
           compressor = "bgzip",
           col_names = !append,
           ...) {
    stopifnot(is.character(file_path) &&
                length(file_path) == 1 &&
                any(endsWith(
                  file_path, c(".gz", ".bed", ".bedGraph")
                )))
    
    args <- list(...)
    
    # User-provided header information
    if (!is.null(header)) {
      header <- paste0("#", header)
    }
    header <-
      paste(c(header, paste0("#", paste(
        colnames(df), collapse = "\t"
      ))),
      collapse = "\n")
    
    if (endsWith(file_path, ".gz")) {
      destfile <- tempfile(fileext = ".bed")
    } else {
      destfile <- file_path
    }
    
    if (col_names)
      cat(paste0(header, "\n"), file = destfile, append = append)
    df %>% write_tsv(
      file = destfile,
      col_names = FALSE,
      # If append is FALSE, bug col_names is TRUE, it means we have to write the header information
      # Thus, we still write in append mode
      append = append || col_names,
      ...
    )
    
    if (endsWith(file_path, ".gz")) {
      system(str_interp("${compressor} ${destfile}"))
      system(str_interp("mv ${destfile}.gz ${file_path}"))
    }
  }