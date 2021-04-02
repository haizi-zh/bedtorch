
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bedtorch

<!-- badges: start -->

[![R-CMD-check](https://github.com/haizi-zh/bedtorch/workflows/R-CMD-check/badge.svg)](https://github.com/haizi-zh/bedtorch/actions)
[![test-coverage](https://github.com/haizi-zh/bedtorch/workflows/test-coverage/badge.svg)](https://github.com/haizi-zh/bedtorch/actions)
<!-- badges: end -->

## Motivation

The goal of bedtorch is to provide a fast and native toolsuite for BED
file manipulation.

In bioinformatics studies, an important type of jobs is related to BED
and BED-like file manipulation. For example, filtering cfDNA fragments
based on length, identifying genomic regions overlapping with certain
genes, filtering features based on epigenetic profiles such as histone
modification levels, finding out SNVs in a certain genomic interval,
etc.

In shell, many of such tasks can be done by some highly-optimized tools,
such as bedtools and bedops. However, sometimes people wants to achieve
this in R environment, not by some external command line tools. This is
especially true when people are developing their methods in R.

R’s native `data.frame` type seems to be a good model for BED-like data,
but in many scenarios it’s peformance is poor, and more importantly, it
lacks many advanced features compared with `dplyr` and `data.table`.
Another choice is using `bedr`, which is bridge linking R codes and
bedtools together. However, this workflow requires lots of disk IO,
since data frames in R needs to be written to temporary folders before
being fed to bedtools. To make things even worse, output from bedtools
are captured and transferred back to R session’s space via `system2`.
This can be extremely slow when the dataset’s size is large. Moreover,
data transferred back to R will lose any information about data type,
index (in the case of `data.table`), column names, etc. Users are on
their own to rebuild the data frame (`tibble` or `data.table`), and
reset data type, column names, etc. This can be very tedious,
error-prone and slow (for rebuilding index). Therefore, it would be
helpful if there is a R package providing these functionalities
natively.

bedtorch is an attempt for this. Under the hood, every BED dataset is
represented as a `data.table` object. Users can apply various operations
on it, such as intersect, merge, etc. Users can also perform any
`data.table` operations as needed.

In terms of the core computation, `bedtorch`’s performance is comparable
to bedtools. However, since no disk IO and data conversion is needed,
users can directly manipulate the dataset in memory, therefore in
practice, via bedtorch, many tasks can be done about one magnitude
faster.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("haizi-zh/bedtorch")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bedtorch)

## Load BED data
bedtbl <-
  read_bed(system.file("extdata", "example2.bed.gz", package = "bedtorch"),
           range = "1:3001-4000")
head(bedtbl)
#>    chrom start  end score1 score2
#> 1:     1  2925 3011    106    181
#> 2:     1  3003 3092     88    193
#> 3:     1  3091 3193    118    212
#> 4:     1  3164 3248     94    211
#> 5:     1  3232 3345    107    205
#> 6:     1  3300 3395     88    193

## Merge intervals, and take the mean of score in each merged group
merged <- merge_bed(bedtbl,
                    operation = list(
                      mean_score = function(x)
                        mean(x$score1)
                    ))
head(merged)
#>    chrom start  end mean_score
#> 1:     1  2925 4016   99.07143

## Find intersections between two datasets
tbl_x <- read_bed(system.file("extdata", "example_merge.bed", package = "bedtorch"))
tbl_y <- read_bed(system.file("extdata", "example_intersect_y.bed", package = "bedtorch"))
head(intersect_bed(tbl_x, tbl_y))
#>    chrom start end score
#> 1:    21    22  25     7
#> 2:    21    26  30     7
#> 3:    21    29  35     9
#> 4:    21    47  49     1
#> 5:    21    47  50     2
#> 6:    21    53  55     5
```

## See also

For more details, refer to the documentation page:
<https://haizi-zh.github.io/bedtorch/>

Many features are inspired by bedtools. Thus, it’s helpful to get
familiar with bedtool’s documentation:
<https://bedtools.readthedocs.io/en/latest/index.html>
