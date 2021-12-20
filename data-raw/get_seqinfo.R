# Prefetch seqinfo for NCBI assemblies

library(tidyverse)

genomes <- data.table::as.data.table(GenomeInfoDb::registered_NCBI_assemblies())
ncbi_seqinfo <- genomes[organism == "Homo sapiens", assembly] %>%
  map(function(assembly) {
    GenomeInfoDb::Seqinfo(genome = assembly)
  })
names(ncbi_seqinfo) <- genomes[organism == "Homo sapiens", assembly]

ncbi_seqinfo$`hs37-1kg` <- ncbi_seqinfo$GRCh37
ncbi_seqinfo$`hs37d5` <- ncbi_seqinfo$GRCh37

usethis::use_data(ncbi_seqinfo, overwrite = TRUE)
