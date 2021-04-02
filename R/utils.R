#' Find a new column name which does not conflict with existing names
#' 
#' @examples
#' # Returns "chrom.1"
#' .available_colname(dt, "chrom")
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