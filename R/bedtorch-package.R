# Home-made essentials tools for bioinformatics

#' @importFrom stringr str_interp str_split str_detect str_sort str_match
#' @importFrom purrr map map_int map_dbl map_lgl map_chr walk pwalk pmap_chr pmap_dfr
#' @importFrom tidyr expand_grid
#' @importFrom dplyr `%>%`
#' @importFrom data.table setnames setkey fread fwrite shift setDT data.table rbindlist as.data.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps mcols `mcols<-` pintersect ranges `ranges<-` width
#' @importFrom S4Vectors from to
NULL

## usethis namespace: start
#' @useDynLib bedtorch, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL