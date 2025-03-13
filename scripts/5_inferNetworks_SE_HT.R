library(NetCoMi)
library(tidyverse)
library(analogue)
library(psych)
library("stringr")

loadGraph <- function(path){
  graph <- read.csv(path, sep=",", header=TRUE)
  return (graph)
}

save.image()
workdir <- snakemake@config$workdir
seed <-snakemake@config$seed
n <-snakemake@config$n
inferthresh <- snakemake@config$inferthresh 

for (sim in snakemake@input[]){
  print(sim)
  tmpOut <- str_replace(str_replace(sim, "abundances","networks"), "filt_abunds.csv", "infmethod.csv")
  #str_replace(tmpOut, "infmethod.csv", "spearman.csv")
  
  M <- loadGraph(sim)
  M['X'] <- NULL
  M <- data.matrix(M, rownames.force = NA)
  
  sim_A <- loadGraph(str_replace(sim, "filt_abunds","filt_base_A"))
  sim_A['X'] <- NULL
  rownames(sim_A) <- colnames(sim_A)
  sim_A <- data.matrix(sim_A, rownames.force = NA)
  
  # SPIECEASI
  net_spieceasi_HT <- netConstruct(M,
                                   dataType = "counts",
                                   measure = "spieceasi",
                                   filtSamp = "none",
                                   filtSampPar = "none",
                                   sparsMethod = 'threshold',
                                   thresh = inferthresh,
                                   weighted = TRUE,
                                   verbose = 3)
  spieceasi_A_HT <- net_spieceasi_HT$adjaMat1
  diag(spieceasi_A_HT) <- 0
  sum(spieceasi_A_HT)
  write.table(spieceasi_A_HT, str_replace(tmpOut, "infmethod.csv", "spieceasi_weighted.csv"), sep=",", dec=".")
  
}