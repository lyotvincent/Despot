# for pipline
setwd("D:/PyCharm/spatial")
h5File <- "h5ads/temp.h5spt"
sptFile <- h5File
name <- "LIBD_151673"
dir <- "data/spatialLIBD/151673"
tenXdir <- dir
ground_truth <- paste0(dir, "/cluster_labels.csv")

pipline_init <-function(sptFile, tenXdir, name = NULL,
                        ground_truth = NULL, ground_sep = ",") {
  Init_spt(sptFile)
  Save_10X_to_spt(tenXdir, sptFile, name = name)
}

pipline_Decontamination <- function(sptFile, tenXdir, method = "SpotClean") {
  if(method == "SpotClean"){
    slide_obj <- Load_10Xh5_to_SpotClean(tenXdir)
    slide_obj <- Decontaminate_SpotClean(slide_obj)
    Save_h5spt_from_SpotClean(slide_obj, sptFile)
  }
}

pipline_Cluster <- function(sptFile, method = "Seurat", h5data = 'matrix'){
  if(method == "BayesSpace"){
    Bayes <- Load_spt_to_SCE(sptFile, h5data = h5data)
    Bayes <- Cluster_BayesSpace(Bayes)
    return(Bayes)
    # Save Bayes Result

  }
  else if(method == "Seurat"){
    seu0 <- Load_spt_to_Seurat(sptFile)
    seu0 <- Preprocess_Seurat(seu0)
    seu0 <- Cluster_Seurat(seu0)
    Save_spt_from_Seurat(sptFile, seu0)
  }
}

pipline_find_hvgs <- function(sptFile, method = "trendsceek"){

}

# Load data to sptFile, and begin the pipline
pipline_init(h5File, dir, name, ground_truth)

# if need, do Decontamination
pipline_Decontamination(h5File, dir)
spt <- h5read(sptFile, '/')
# Do cluster
Bayes <- pipline_Cluster(h5File, method = "BayesSpace", h5data = "spotClean_mat")

# find hvgs


# load in scRNA-seq Data, and do Deconvolution

# if need, do Enhance

# Do plot, Save all the data into spt
