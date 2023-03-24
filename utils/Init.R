# Title     : Initialization of sptFile
# Objective : Pip
# Created by: rzh
# Created on: 2022/3/10

source("sptranr/R/transpar.R")

# decoding params
params <- fromJSON(file = "params.json")

# set working directory
setwd(params$working_dir)

# set data path
dataPath <- params$dataPath

# set sptFile
sptFile <- params$sptFile

# get platform
platform <- params$platform

# set name
name <- params$name

# has ground_truth?
ground.truth <- paste0(params$working_dir, '/', dataPath, '/', params$ground_truth)

ground.name <- params$ground_name

# Init sptFile
Init_spt(sptFile)

# Save important configs to sptFile
Save_cfg_to_spt(sptFile, params)

if(platform %in% c("10X", "10X_Visium", "Visium", "visium", "10X Visium")){
  # is filtered_matrix?
  filtered.matrix <- params$filter_matrix

  # save hires.png?
  save.hires <- params$load_hires

  Save_10X_to_spt(dir = dataPath,
                  sptFile = sptFile,
                  filtered.matrix = filtered.matrix,
                  name = name,
                  ground.truth = ground.truth,
                  ground.name = ground.name,
                  save.hires = save.hires)
} else if(platform %in% c("Slide-seq", "Slide-seqV2", "slide-seq", 'slide-seqV2')) {
  Save_SlideSeq_to_spt(dir = dataPath,
                 sptFile = sptFile,
                 name = name)
} else if(platform %in% c("ST", "st")){
  dataName <- params$dataName
  st <- paste0(dataPath, '/', dataName)
  Save_ST_to_spt(st = st,
                 sptFile = sptFile,
                 name = name)
}

