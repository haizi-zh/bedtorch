# Prefetch seqinfo for Homo sapiens assemblies

library(tidyverse)

genomes <- rbind(
  data.table::as.data.table(GenomeInfoDb::registered_NCBI_assemblies())[organism == "Homo sapiens", .(organism, assembly)],
  data.table::as.data.table(GenomeInfoDb::registered_UCSC_genomes())[organism == "Homo sapiens", .(organism, assembly = genome)]
)
hs_seqinfo <- genomes[, assembly] %>%
  map(function(assembly) {
    cat(str_interp("Loading assembly ${assembly} ... \n"))
    GenomeInfoDb::Seqinfo(genome = assembly)
  })
names(hs_seqinfo) <- genomes[, assembly]
hs_seqinfo$`hs37-1kg` <- hs_seqinfo$GRCh37
hs_seqinfo$`hs37d5` <- hs_seqinfo$GRCh37

usethis::use_data(hs_seqinfo, overwrite = TRUE)
