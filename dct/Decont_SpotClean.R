# Title     : Using SpotClean for Decontamination
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_SpotClean.R")

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
dataPath <- params$dataPath
platform <- params$platform
if(platform %in% c("10X", "10X_Visium", "Visium", "visium", "10X Visium")){
  slide_obj <- Load_10Xh5_to_SpotClean(dataPath, smdFile)
  slide_obj <- Decontaminate_SpotClean(slide_obj)
  Save_smd_from_SpotClean(smdFile, slide_obj)
} else{
  message(paste0("SpotClean only support 10X series. But your platform is ", platform))
}


