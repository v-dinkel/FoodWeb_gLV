# Some installations which might need to be handled manually:
#install.packages("remotes")
#install.packages("ragg")
#install.packages("devtools")
#install.packages("BiocManager")
#remotes::install_version("Matrix", version = "1.6")
#install.packages("pulsar")
#devtools::install_github("zdk123/SpiecEasi")
#install.packages("sn")
#devtools::install_github("GraceYoon/SPRING")
#remotes::install_version("Hmisc", version = "5.1")
#devtools::install_github("stefpeschel/NetComi", repos= c("https://cloud.r-project.org/", BiocManager::repositories()))
#install.packages("RcppGSL") #this is system requirement
#install.packages("ecoCopula")
#devtools::install_github("tpq/propr")

library(NetCoMi)
library(tidyverse)
library(analogue)
library(psych)
library(ecoCopula)
library("stringr")

# helper function to load a csv
loadGraph <- function(path){
  graph <- read.csv(path, sep=",", header=TRUE)
  return (graph)
}

save.image()

#read snakemake config variables 
workdir <- snakemake@config$workdir
seed <-snakemake@config$seed
n <-snakemake@config$n
inferthresh <- snakemake@config$inferthresh 

# iterate over input files from snakemake
for (sim in snakemake@input[]){
  print(sim)
  tmpOut <- str_replace(str_replace(sim, "abundances","networks"), "filt_abunds.csv", "infmethod.csv")
  
  # load the 
  M <- loadGraph(sim)
  M['X'] <- NULL
  M <- data.matrix(M, rownames.force = NA)
  
  sim_A <- loadGraph(str_replace(sim, "filt_abunds","filt_base_A"))
  sim_A['X'] <- NULL
  rownames(sim_A) <- colnames(sim_A)
  sim_A <- data.matrix(sim_A, rownames.force = NA)
  
  # set the threshold to default 0.4. This is done in all methods, but may be disregarded for methods which have other sparsification methods.
  thresh_spearman <- inferthresh
  thresh_propr <- inferthresh
  thresh_spieceasi <- inferthresh
  thresh_ccrepe <- inferthresh
  thresh_sparcc <- inferthresh
  
  # SPEARMAN
  net_spearman <- netConstruct(M,
                               dataType = "counts",
                               measure = "spearman",
                               filtSamp = "none",
                               filtSampPar = "none",
                               normMethod = "TSS",
                               weighted = FALSE,
                               sparsMethod = 'threshold',
                               thresh = thresh_spearman,
                               verbose = 3)
  spearman_A <- net_spearman$adjaMat1
  diag(spearman_A) <- 0
  write.table(spearman_A, str_replace(tmpOut, "infmethod.csv", "spearman.csv"), sep=",", dec=".")
  
  # PROPR
  net_propr <- netConstruct(M,
                            dataType = "counts",
                            measure = "propr",
                            filtSamp = "none",
                            filtSampPar = "none",
                            sampleSize = 100,
                            sparsMethod = 'threshold',
                            normMethod = "TSS",
                            weighted = FALSE,
                            thresh = thresh_propr,
                            verbose = 3)
  propr_A <- net_propr$adjaMat1
  diag(propr_A) <- 0
  write.table(propr_A, str_replace(tmpOut, "infmethod.csv", "propr.csv"), sep=",", dec=".")
  
  # SPIECEASI
  net_spieceasi <- netConstruct(M,
                                dataType = "counts",
                                measure = "spieceasi",
                                filtSamp = "none",
                                filtSampPar = "none",
                                normMethod = "TSS",
                                sparsMethod = 'threshold',
                                weighted = FALSE,
                                thresh = thresh_spieceasi,
                                verbose = 3)
  spieceasi_A <- net_spieceasi$adjaMat1
  diag(spieceasi_A) <- 0
  write.table(spieceasi_A, str_replace(tmpOut, "infmethod.csv", "spieceasi.csv"), sep=",", dec=".")
  
  # CCREPE
  net_ccrepe <- netConstruct(M,
                             dataType = "counts",
                             measure = "ccrepe",
                             filtSamp = "none",
                             filtSampPar = "none",
                             normMethod = "TSS",
                             sparsMethod = 'threshold',
                             sampleSize = 100,
                             weighted = FALSE,
                             thresh = thresh_ccrepe,
                             verbose = 3)
  ccrepe_A <- net_ccrepe$adjaMat1
  diag(ccrepe_A) <- 0
  write.table(ccrepe_A, str_replace(tmpOut, "infmethod.csv", "ccrepe.csv"), sep=",", dec=".")
  
  # SPARCC
  net_sparcc <- netConstruct(M,
                             dataType = "counts",
                             measure = "sparcc",
                             filtSamp = "none",
                             filtSampPar = "none",
                             sparsMethod = 'threshold',
                             normMethod = "TSS",
                             weighted = FALSE,
                             thresh = thresh_sparcc,
                             verbose = 3)
  sparcc_A <- net_sparcc$adjaMat1
  diag(sparcc_A) <- 0
  write.table(sparcc_A, str_replace(tmpOut, "infmethod.csv", "sparcc.csv"), sep=",", dec=".")
  
  # ECOCOPULA
  M_rnd <- ceiling( M )
  tmp <- as.data.frame(rep(1:1, each=length(M_rnd[,1])))
  colnames(tmp) <- "tmp"
  M_rnd_mod <- stackedsdm(M_rnd, ~., data = tmp, ncores=4)
  M_ecocopula = cgr( M_rnd_mod, lambda = 0.5 )
  M_ecocopula_A = M_ecocopula$best_graph$graph
  diag(M_ecocopula_A) <- 0
  write.table(M_ecocopula_A, str_replace(tmpOut, "infmethod.csv", "ecocopula.csv"), sep=",", dec=".")
}
