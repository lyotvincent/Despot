# Title     : generate Normal scRNA-seq Data
# Objective : Pip
# Created by: rzh
# Created on: 2022/4/25

source("sptranr/R/transpar.R")
source("sptranr/R/_Seurat.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
dataPath <- params$dataPath
scdir <- params$scDir
scfiles <- params$scFiles

seu <- Load_txtsc_to_Seurat(scfiles, scdir)

seu <- Analysis_scRNA_seq_Seurat(seu, resolution = 0.8)

sce <- as.SingleCellExperiment(seu)
rm(seu)
sce <- Analysis_scRNA_seq(sce)
Save_scRNAseq_to_spt(sptFile, sce)
