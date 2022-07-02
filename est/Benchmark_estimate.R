# Title     : Benchmark for spTRS
# Objective : Pip
# Created by: rzh
# Created on: 2022/4/8

source("sptranr/R/transpar.R")
source("sptranr/R/_benchmark.R")
source("sptranr/R/_Seurat.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
imgdir <- paste0(params$dataPath, "/spatial")

mark <- Benchmark_estimate(sptFile, h5data)
Save_spt_from_Benchmark_est(sptFile, h5data, mark)
