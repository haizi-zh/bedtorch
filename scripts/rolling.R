devtools::load_all()
library(tidyverse)

# Performance ----

## rollsum ----

n <- 50e6L
k <- 500
test_data <- runif(n)
system.time(zoo::rollsum(test_data, k))
system.time(bedtorch::rollsum(test_data, k))

system.time(zoo::rollmean(test_data, k))
system.time(bedtorch::rollmean(test_data, k))
