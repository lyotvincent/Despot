# Title     : Using BayesSpace for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_BayesSpace.R")
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

for(decont in params$Decontamination){
  h5data <- Create_spt_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }
  Bayes <- Load_spt_to_SCE(sptFile, h5data = h5data)
  Bayes <- Cluster_BayesSpace(Bayes)
  Save_spt_from_BayesSpace(sptFile, Bayes, save.enhanced = F, h5data = h5data)
}


