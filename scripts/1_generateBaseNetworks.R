save.image()
workdir <- snakemake@config$workdir
seed <-snakemake@config$seed
n <-snakemake@config$n
nettype <-snakemake@config$nettype
options(warn=-1)

library(devtools)
library(SpiecEasi)

generateSyntheticNetworks <- function(outpath, type) {
  set.seed(seed)
  graph <- make_graph(type, n, n)
  write.csv(graph, outpath, row.names = TRUE)
}

generateSyntheticNetworks(snakemake@output[[1]], nettype)