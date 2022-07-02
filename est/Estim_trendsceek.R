# Title     : Using trendsceek for estimation
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/12

source("sptranr/R/transpar.R")
source("sptranr/R/_trendsceek.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
h5data <- "matrix"
if(params$Decontamination == "SpotClean"){
  h5data <- "SpotClean_mat"
}
if(params$Decontamination == "SPCS"){
  h5data <- "SPCS_mat"
}
imgdir <- paste0(params$dataPath, "/spatial")

seu <- Load_spt_to_Seurat(sptFile,imgdir, h5data)
seu <- Preprocess_Seurat(seu)
