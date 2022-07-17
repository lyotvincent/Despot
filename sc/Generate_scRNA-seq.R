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

# Load scRNA-seq data format by users' definition
scbar <- paste0(scdir, '/', params$scBarcodes)
scfea <- paste0(scdir, '/', params$scFeatures)
scmat <- paste0(scdir, '/', params$scMatrix)
scgth <- paste0(scdir, '/', params$scGroundTruth)
scgnm <- params$scGroundName

scType <- params$scType
if(scType == "Any"){
  sce <- Load_sc_to_SCE(scbar, scfea, scmat, scgth, scgnm)
  sce <- Analysis_scRNA_seq(sce)
} else  if(scType == "10X"){
  sce <- Load_10Xsc_to_SCE(scdir)
  sce <- Analysis_scRNA_seq_noanno(sce, anno.signature = scMarkers)
}

Save_scRNAseq_to_spt(sptFile, sce)
