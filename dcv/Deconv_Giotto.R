source("sptranr/R/transpar.R")
source("sptranr/R/_Giotto.R")
source("sptranr/R/_scRNA-seq.R")
library(tidyverse)
# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
python_path <- params$pythonPath
temp_dir <- params$tempdir

sce <- Load_sptsc_to_SCE(sptFile, h5data = "scRNA_seq")
sce <- Analysis_scRNA_seq(sce)
sce0 <- sce[sce@metadata$HVGs, ]
Sig <- aggregate(data.frame(t(as.matrix(sce0@assays@data$logcounts))), by=list(sce0$free_annotation), FUN=mean)
rownames(Sig) <- Sig[[1]]
Sig <- Sig[, -1]
Sig0 <- t(as.matrix(Sig))

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
  gobj <- runDWLSDeconv(gobj, sign_matrix = Sig0)
  Save_spt_from_Giotto_dcv(sptFile, gobj, h5data)
}
