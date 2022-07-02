# Title     : Using SPCS for Decontamination
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")
source("sptranr/R/_SPCS.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile

sce <- Load_spt_to_SCE(sptFile)
SPCS_obj <- Decontaminate_SPCS(sce)
Save_spt_from_SPCS(sptFile, SPCS_obj)
