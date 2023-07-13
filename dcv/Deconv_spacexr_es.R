# Title     : Using spacexr for Deconvolution
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/18

source("sptranr/R/transpar.R")
source("sptranr/R/_spacexr.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
sptFile <- params$sptFile
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
  # read References
  sce <- Load_sptsc_to_SCE(sptFile, "sc-ref-es")
  anno <- sce$free_annotation
  sce$free_annotation <- gsub("/", "^", anno)   #cell-types in RCTD don't support `/`, transport to `^`
  anno <- table(anno)
  sr <- Load_spt_to_SpatialRNA(sptFile, h5data)
  ref <- GenerateRef_spacexr(sce)
  rctd <- Deconvolution_spacexr(sr, ref)
  rctd@results$weights@Dimnames[[2]] <- attr(anno, "names")
  Save_spt_from_spacexr(sptFile, h5data, rctd, pp_mtd = "es")
}
