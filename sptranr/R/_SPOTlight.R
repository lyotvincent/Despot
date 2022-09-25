library(SPOTlight)

# Load SCE data
Deconvolution_SPOTlight <- function(spe, sce){

  res <- SPOTlight(
    x = sce,
    y = spe,
    assay = "counts",
    groups = sce$free_annotation,
    hvg = sce@metadata$HVGs,
    mgs = sce@metadata$marker_genes,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")

  return(res)
}

# save spt from SPOTlight results
Save_spt_from_SPOTlight <- function(sptFile, h5data, spotlight, pp_mtd="raw"){

  results <- as(spotlight$mat, "dgeMatrix")

  if(pp_mtd == "raw"){
    mtd_name <- "SPOTlight"
  } else if(pp_mtd == "es"){
    mtd_name <- "SPOTlight_es"
  } else if(pp_mtd == "vae"){
    mtd_name <- "SPOTlight_vae"
  }

  h5createGroup(sptFile, paste0(h5data, '/deconv/', mtd_name))

  # save 1d weights
  Create_spt_array1d(sptFile,
                     arr = results@x,
                     sptloc = paste0(h5data, '/deconv/', mtd_name, '/weights'),
                     mode = 'double')
  # save shape
  Create_spt_array1d(sptFile,
                     arr = results@Dim,
                     sptloc = paste0(h5data, '/deconv/', mtd_name ,'/shape'),
                     mode = 'integer')
  # save dim names
  Create_spt_array1d(sptFile,
                     arr = rownames(results),
                     sptloc = paste0(h5data, '/deconv/', mtd_name ,'/barcodes'),
                     mode = 'character')
  Create_spt_array1d(sptFile,
                     arr = colnames(results),
                     sptloc = paste0(h5data, '/deconv/', mtd_name ,'/cell_type'),
                     mode = 'character')
}
