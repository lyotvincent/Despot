# Title     : Using SpotClean for Decontamination
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_SpotClean.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
dataPath <- params$dataPath

slide_obj <- Load_10Xh5_to_SpotClean(dataPath, sptFile)
slide_obj <- Decontaminate_SpotClean(slide_obj)
Save_spt_from_SpotClean(sptFile, slide_obj)

