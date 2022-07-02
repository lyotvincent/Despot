# Title     : Using spacexr for Deconvolution
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/18

source("sptranr/R/transpar.R")
source("sptranr/R/_spacexr.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
refdir <- paste0(params$tempdir, "/references")

for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  # read References
  sce <- Load_sptsc_to_SCE(sptFile, "scRNA_seq")
  sr <- Load_spt_to_SpatialRNA(sptFile)
  ref <- GenerateRef_spacexr(sce)
  rctd <- Deconvolution_spacexr(sr, ref)
  Save_spt_from_spacexr(sptFile, h5data, rctd)
}
