# Title     : generate Normal scRNA-seq Data
# Objective : Pip
# Created by: rzh
# Created on: 2022/4/25

source("sptranr/R/transpar.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
dataPath <- params$dataPath
scdir <- params$scDataPath
scMarkers <- params$scMarkers

sce <- Load_10Xsc_to_SCE(scdir)

sce <- Analysis_scRNA_seq_noanno(sce, anno.signature = scMarkers)

Save_scRNAseq_to_spt(sptFile, sce)
