# Title     : Using scater to generate TabulaData
# Objective : Pip
# Created by: rzh
# Created on: 2022/4/10

source("sptranr/R/transpar.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
dataPath <- params$dataPath

# read TabulaData
TabulaType <- params$TabulaType
TabulaData <- params$TabulaData

sce <- Load_TabulaData_to_SCE(TabulaType, TabulaData)

sce <- Analysis_TabulaData(sce)

Save_scRNAseq_to_spt(sptFile, sce)
