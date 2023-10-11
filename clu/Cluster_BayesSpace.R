# Title     : Using BayesSpace for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_BayesSpace.R")
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

for(decont in params$Decontamination){
  h5data <- Create_smd_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }
  Bayes <- Load_smd_to_SCE(smdFile, h5data = h5data)
  Bayes <- Cluster_BayesSpace(Bayes)
  Save_smd_from_BayesSpace(smdFile, Bayes, save.enhanced = F, h5data = h5data)
}


