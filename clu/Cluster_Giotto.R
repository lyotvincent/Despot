# Title     : Using gobjrat for clustering
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/25

source("sptranr/R/transpar.R")
source("sptranr/R/_Giotto.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
python_path <- params$pythonPath
temp_dir <- params$tempdir

for(decont in params$Decontamination){
  h5data <- "matrix"
  if(decont == "SpotClean"){
    h5data <- "SpotClean_mat"
  }
  if(decont == "SPCS"){
    h5data <- "SPCS_mat"
  }

  gobj <- Load_spt_to_Giotto(sptFile, h5data, python_path, temp_dir)
  gobj <- Preprocess_Giotto(gobj)
  gobj <- Cluster_Giotto(gobj)
  Save_spt_from_Giotto_clu(sptFile, gobj, h5data)
}
