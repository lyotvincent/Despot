# checking whether SpotClean is in the library
source("sptranr/R/_Loading.R")

Check_Load_GithubPackages("SpotClean", "zijianni/SpotClean@4fd969c")
Check_Load_BiocPackages("rhdf5")

Load_10Xh5_to_SpotClean <- function(dataPath, smdFile){
  raw_path_h5 <- paste0(dataPath, "/raw_feature_bc_matrix.h5")
  raw_path <- paste0(dataPath, "/raw_feature_bc_matrix")
  csv_path <- paste0(dataPath, "/spatial/tissue_positions_list.csv")
  img_path <- paste0(dataPath, "/spatial/tissue_lowres_image.png")
  fac_path <- paste0(dataPath, "/spatial/scalefactors_json.json")
  # read feature matrix
  if(file.exists(raw_path_h5)){
    raw_mat <- Read10xRawH5(h5_file = raw_path_h5)
  }
  else if(file.exists(raw_path)){
    raw_mat <- Read10xRaw(count_dir = raw_path)
  }
  else{
    message("SpotClean needs raw file in 10X Visium data.")
    return(0)
  }
  # read slide info
  slide_info <- Read10xSlide(tissue_csv_file = csv_path,
                             tissue_img_file = img_path,
                             scale_factor_file = fac_path)
  # load filterd barcodes
  filtered_barcodes <- rhdf5::h5read(smdFile, 'matrix/barcodes')
  in_tissue <- as.integer(slide_info$slide$barcode %in% filtered_barcodes)
  slide_info$slide$tissue <- in_tissue
  # Create slide object
  slide_obj <- CreateSlide(raw_mat, slide_info)
  return(slide_obj)
}

Decontaminate_SpotClean <- function(slide_obj,
                                    maxit = 10,
                                    candidate_radius = 20){
  # Decontaminate the data
  slide_obj <- SpotClean(slide_obj, maxit=10, candidate_radius = 20)
  return(slide_obj)
}

Save_smd_from_SpotClean <- function(smdFile, slide_obj){
  # the ARC score is a conserved lower bound of the proportion
  # of contamination in observed tissue spots
  arc <- slide_obj@metadata$ARC_score

  h5data <- "SpotClean_mat"
  rhdf5::h5createGroup(smdFile, h5data)
  Create_smd_array1d(smdFile,
                     arr = arc,
                     smdloc = paste0(h5data, "/ARCscore"),
                     mode = "double")

  slide_obj <- slide_obj[,order(colnames(slide_obj))]

  # make the number of barcods equal to the filtered matrix
  barcodes <- rhdf5::h5read(smdFile, 'matrix/barcodes')
  slide_obj <- slide_obj[, colnames(slide_obj) %in% barcodes]
  # create barcodes
  Create_smd_array1d(smdFile,
                     arr = slide_obj@colData@rownames,
                     smdloc = paste0(h5data, "/barcodes"),
                     mode = "character")

  # create data, indices, indptr, shape
  Create_smd_array1d(smdFile,
                     arr = slide_obj@assays@data@listData$decont@x,
                     smdloc = paste0(h5data, "/data"),
                     mode = "double")
  Create_smd_array1d(smdFile,
                     arr = slide_obj@assays@data@listData$decont@i,
                     smdloc = paste0(h5data, "/indices"),
                     mode = "integer")
  Create_smd_array1d(smdFile,
                     arr = slide_obj@assays@data@listData$decont@p,
                     smdloc = paste0(h5data, "/indptr"),
                     mode = "integer")
  Create_smd_array1d(smdFile,
                     arr = slide_obj@assays@data@listData$decont@Dim,
                     smdloc = paste0(h5data, "/shape"),
                     mode = "integer")

  # create features
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/features"))
  Create_smd_array1d(smdFile,
                     arr = slide_obj@NAMES,
                     smdloc = paste0(h5data, "/features/name"),
                     mode = "character")
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/features/is_HVG"))

  # create contamination rate
  Create_smd_array1d(smdFile,
                     arr = slide_obj@metadata$contamination_rate,
                     smdloc = paste0(h5data, "/contamination_rate"),
                     mode = "double")

  # create `idents`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/idents"))
  # create `deconv`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/deconv"))
  # create `benchmark`
  rhdf5::h5createGroup(smdFile, paste0(h5data, "/benchmark"))
}
