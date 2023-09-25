library(CARD)
library(MuSiC)
source("sptranr/R/transpar.R")
source("sptranr/R/_scRNA-seq.R")

# install TOAST, MuSiC


Deconvolution_CARD <- function(spe, sce){

  rownames(sce@assays@data$counts) <- Make_unique(rownames(sce))
  colnames(sce@assays@data$counts) <- Make_unique(colnames(sce))
  sc_count <- sce@assays@data$counts
  sc_meta <- data.frame(cellType=sce$free_annotation,
                        sampleInfo = 'sample1',
                        row.names = colnames(sc_count))
  spatial_location <- data.frame(x=spe@colData[['row']],
                            y=spe@colData[['col']],
                            row.names = colnames(spe))
  spatial_count <- spe@assays@data$counts
  card <- createCARDObject(sc_count=sc_count,
                           sc_meta = sc_meta,
                           spatial_count = spatial_count,
                           spatial_location = spatial_location,
                           ct.varname = 'cellType',
                           ct.select = unique(sc_meta$cellType),
                           sample.varname = "sampleInfo",
                           minCountSpot = 5,
                           minCountGene = 0)

  card <- CARD_deconvolution(card)
  return(card@Proportion_CARD)
}

# save smd from CARD results
Save_smd_from_CARD <- function(smdFile, h5data, card, pp_mtd="raw"){

  results <- as(card, "dMatrix")

  if(pp_mtd == "raw"){
    mtd_name <- "CARD"
  } else if(pp_mtd == "es"){
    mtd_name <- "CARD_es"
  } else if(pp_mtd == "vae"){
    mtd_name <- "CARD_vae"
  }

  h5createGroup(smdFile, paste0(h5data, '/deconv/', mtd_name))

  # save 1d weights
  Create_smd_array1d(smdFile,
                     arr = results@x,
                     smdloc = paste0(h5data, '/deconv/', mtd_name, '/weights'),
                     mode = 'double')
  # save shape
  Create_smd_array1d(smdFile,
                     arr = results@Dim,
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/shape'),
                     mode = 'integer')
  # save dim names
  Create_smd_array1d(smdFile,
                     arr = rownames(results),
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/barcodes'),
                     mode = 'character')
  Create_smd_array1d(smdFile,
                     arr = colnames(results),
                     smdloc = paste0(h5data, '/deconv/', mtd_name ,'/cell_type'),
                     mode = 'character')
}
