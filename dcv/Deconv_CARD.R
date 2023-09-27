source("sptranr/R/transpar.R")
source("sptranr/R/_CARD.R")
source("sptranr/R/_scRNA-seq.R")

# decoding params
params <- fromJSON(file = "params.json")
smdFile <- params$smdFile
platform <- params$platform
if(is.null(platform)){
  platform <- "10X_Visium"  # default using 10X_Visium
}
sce <- Load_smdsc_to_SCE(smdFile, "scRNA_seq")
for(decont in params$Decontamination){
  h5data <- Create_smd_h5data(decont)
  if(platform != "10X_Visium" && h5data == "SpotClean_mat"){
    message("SpotClean only Support 10X Visium data, skip it.")
    next
  }

  # read References
  spe <- Load_smd_to_SCE(smdFile, h5data)
  card <- Deconvolution_CARD(spe, sce)
  Save_smd_from_CARD(smdFile, h5data, card, pp_mtd = "raw")
}
