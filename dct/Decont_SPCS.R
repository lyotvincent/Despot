# Title     : Using SPCS for Decontamination
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_SPCS.R")

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile

sce <- Load_smd_to_SCE(smdFile)
SPCS_obj <- Decontaminate_SPCS(sce)
Save_smd_from_SPCS(smdFile, SPCS_obj)
