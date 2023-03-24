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
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}

for(decont in params$Decontamination){
  h5data <- Create_spt_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }
  gobj <- Load_spt_to_Giotto(sptFile, h5data, python_path, temp_dir)
  gobj <- Preprocess_Giotto(gobj)
  gobj <- Cluster_Giotto(gobj)
  Save_spt_from_Giotto_clu(sptFile, gobj, h5data)
}
