# Title     : Using BayesSpace for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_BayesSpace.R")
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }
  Bayes <- Load_spt_to_SCE(sptFile, h5data = h5data)
  Bayes <- Cluster_BayesSpace(Bayes)
  Save_spt_from_BayesSpace(sptFile, Bayes, save.enhanced = F, h5data = h5data)
}


