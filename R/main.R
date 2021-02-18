# Home-made essentials tools for bioinformatics

#' @import readr
#' @import stringr
#' @import dplyr


# This function loads a genbed file, which is simply a BED-like tsv file,
# with a header line (optional) containing column names, and the types of
# the first three columns are guaranteed to be "cii"
#
# col_names: can be NULL (infer from header), or a character vector
# (user-provided column names).
#
#' @export
load_genbed <-
  function(file_path,
           col_names = NULL,
           na = NULL,
           col_types = NULL,
           ...) {
    if (any(startsWith(file_path, c("http://", "https://"))) &
        endsWith(file_path, ".gz")) {
      tmpbed <- tempfile(fileext = ".gz")
      download.file(url = file_path, destfile = tmpbed)
      genbed <-
        load_genbed(
          tmpbed,
          col_names = col_names,
          na = na,
          col_types = col_types,
          ...
        )
      unlink(tmpbed)
      return(genbed)
    }
    
    # Try identify column names from the header
    get_col_names <- function() {
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
    
    if (is.null(col_names)) {
      col_names <- get_col_names()
    }
    
    if (is.null(col_types)) {
      col_types <-
        paste(c("cii", rep_len("?", length(col_names) - 3)), collapse = "")
    }
    
    # string for missing values
    na <- c(c("", "NA"), na) %>% unique()
    
    genbed <- read_tsv(
      file = file_path,
      comment = "#",
      col_names = col_names,
      col_types = col_types,
      na = na
    )
    
    return(genbed)
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