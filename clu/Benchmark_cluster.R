# Title     : Benchmark for spTRS
# Objective : Pip
# Created by: rzh
# Created on: 2022/4/8

source("sptranr/R/transpar.R")
source("sptranr/R/_benchmark.R")

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
imgdir <- paste0(params$dataPath, "/spatial")
for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  mark <- Benchmark_Clustering(smdFile, h5data)
  Save_smd_from_Benchmark_clu(smdFile, h5data, mark)
  message(paste0("Benchmark for clustering has been saved in ", h5data, "/benchmark/cluster."))
}
