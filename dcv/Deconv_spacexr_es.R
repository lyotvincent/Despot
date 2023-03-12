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
params <- h5read(sptFile, "configs")

for(decont in params$Decontamination){
  h5data <- Create_spt_h5data(decont)

  # read References
  sce <- Load_sptsc_to_SCE(sptFile, "sc-ref-es")
  sr <- Load_spt_to_SpatialRNA(sptFile, h5data)
  ref <- GenerateRef_spacexr(sce)
  rctd <- Deconvolution_spacexr(sr, ref)
  Save_spt_from_spacexr(sptFile, h5data, rctd, pp_mtd = "es")
}
