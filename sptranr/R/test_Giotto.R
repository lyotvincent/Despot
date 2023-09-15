library(dplyr)
library(Giotto)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(Seurat)

temp_dir = 'mini_visium/'
myinstructions = createGiottoInstructions(save_dir = temp_dir,
                                          save_plot = FALSE,
                                          show_plot = F,
                                          python_path = '~/anaconda3/envs/spatial/bin/python')

expr_path = system.file("extdata", "visium_DG_expr.txt.gz", package = 'Giotto')
loc_path = system.file("extdata", "visium_DG_locs.txt", package = 'Giotto')
mini_visium <- createGiottoObject(raw_exprs = expr_path,
                                  spatial_locs = loc_path,
                                  instructions = myinstructions)

## 1. read image
png_path = system.file("extdata", "deg_image.png", package = 'Giotto')
mg_img = magick::image_read(png_path)
mg_img
## 2. test and modify image alignment
mypl = spatPlot(mini_visium, return_plot = T, point_alpha = 0.8)
orig_png = createGiottoImage(gobject = mini_visium, mg_object = mg_img, name = 'image',
                             xmax_adj = 450, xmin_adj = 550,
                             ymax_adj = 200, ymin_adj = 200)
mypl_image = addGiottoImageToSpatPlot(mypl, orig_png)
mypl_image
## 3. add images to Giotto object ##
image_list = list(orig_png)
mini_visium = addGiottoImage(gobject = mini_visium,
                             images = image_list)
showGiottoImageNames(mini_visium)

# explore gene and cell distribution
p1 <- filterDistributions(mini_visium, detection = 'genes') + ggtitle(" Gene  Distributions")
p2 <-  filterDistributions(mini_visium, detection = 'cells')+ ggtitle(" Cell  Distributions")
p1 + p2

filterCombinations(mini_visium,
                   expression_thresholds = c(1),
                   gene_det_in_min_cells = c(20, 20, 50, 50),
                   min_det_genes_per_cell = c(100, 200, 100, 200))

# filter and normalize
mini_visium <- filterGiotto(gobject = mini_visium,
                            expression_threshold = 1,
                            gene_det_in_min_cells = 50,
                            min_det_genes_per_cell = 100,
                            expression_values = c('raw'),
                            verbose = T)
mini_visium <- normalizeGiotto(gobject = mini_visium, scalefactor = 6000, verbose = T)
mini_visium <- addStatistics(gobject = mini_visium)


mini_visium <- calculateHVG(gobject = mini_visium)

# runPCA
mini_visium <- runPCA(gobject = mini_visium)
screePlot(mini_visium, ncp = 30)

ppca <- plotPCA(gobject = mini_visium)
mini_visium <- runUMAP(mini_visium, dimensions_to_use = 1:10)
pumap<- plotUMAP(gobject = mini_visium)
mini_visium <- runtSNE(mini_visium, dimensions_to_use = 1:10)
ptsne <- plotTSNE(gobject = mini_visium)

mini_visium <- createNearestNetwork(gobject = mini_visium, dimensions_to_use = 1:10, k = 15)
mini_visium <- doLeidenCluster(gobject = mini_visium, resolution = 0.4, n_iterations = 1000)

# visualize UMAP cluster results
plotUMAP(gobject = mini_visium, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)

# visualize UMAP and spatial results
spatDimPlot(gobject = mini_visium, cell_color = 'leiden_clus', spat_point_shape = 'voronoi')

# spatial voronoi plot with selected clusters
p1 <- spatPlot(mini_visium, point_shape = 'voronoi', cell_color ='leiden_clus', select_cell_groups = c(1,2,3))

# spatial voronoi plot without showing not selected clusters
p2 <- spatPlot(mini_visium, point_shape = 'voronoi', cell_color ='leiden_clus', select_cell_groups = c(1,2,3), show_other_cells = F)

# spatial voronoi plot without showing not selected cells, but showing the voronoi borders
p3 <- spatPlot(mini_visium, point_shape = 'voronoi', cell_color ='leiden_clus', select_cell_groups = c(1,2,3,4), show_other_cells = F, vor_border_color = 'black')



p1+p2+ p3

mast_markers = findMarkers_one_vs_all(gobject = mini_visium,
                                      method = 'mast',
                                      expression_values = 'normalized',
                                      cluster_column = 'leiden_clus')


topgenes_mast <- mast_markers[, head(.SD, 2), by = 'cluster']$genes  # 这里的.SD 用的很溜啊
topgenes_mast
violinPlot(mini_visium, genes = topgenes_mast, cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right')

topgenes_mast = mast_markers[, head(.SD, 6), by = 'cluster']$genes
topgenes_mast

plotMetaDataHeatmap(mini_visium, selected_genes = topgenes_mast,
                    metadata_cols = c('leiden_clus'))

clusters_cell_types = paste0(mast_markers[, head(.SD, 1), by = 'cluster']$genes,"_cells")
clusters_cell_types
names(clusters_cell_types) = 1:length(unique(mini_visium@cell_metadata$leiden_clus))
mini_visium@cell_metadata
mini_visium = annotateGiotto(gobject = mini_visium, annotation_vector = clusters_cell_types,
                             cluster_column = 'leiden_clus', name = 'cell_types')

# check new cell metadata
pDataDT(mini_visium)

## cell type signatures ##
## combination of all marker genes identified
sign_matrix_path = system.file("extdata", "sig_matrix.txt", package = 'Giotto')
brain_sc_markers = data.table::fread(sign_matrix_path) # file don't exist in data folder
sig_matrix = as.matrix(brain_sc_markers[,-1]); rownames(sig_matrix) = brain_sc_markers$Event
sig_matrix[1:4,1:4]

mini_visium = runSpatialEnrich(mini_visium,
                               sign_matrix = sig_matrix,
                               enrich_method = 'PAGE') #default = 'PAGE'

## heatmap of enrichment versus annotation (e.g. clustering result)
cell_types = colnames(sig_matrix)
plotMetaDataCellsHeatmap(gobject = mini_visium,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8, y_text_size = 8)

enrichment_results = mini_visium@spatial_enrichment$PAGE
enrich_cell_types = colnames(enrichment_results)
enrich_cell_types = enrich_cell_types[enrich_cell_types != 'cell_ID']

enrich_cell_types
## spatplot
?spatCellPlot
spatCellPlot(gobject = mini_visium, spat_enr_names = 'PAGE',
             cell_annotation_values = enrich_cell_types,show_image=T,show_legend=F,
             cow_n_col = 3,coord_fix_ratio = NULL, point_size = 1)

mini_visium <- createSpatialGrid(gobject = mini_visium,
                                 sdimx_stepsize = 300,
                                 sdimy_stepsize = 300,
                                 minimum_padding = 50)
showGrids(mini_visium)
?spatPlot
spatPlot(gobject = mini_visium, show_grid = T, point_size = 1.5,
         show_image=T,group_by="cell_types")

annotated_grid = annotateSpatialGrid(mini_visium)
annotated_grid

annotated_grid_metadata = annotateSpatialGrid(mini_visium,
                                              cluster_columns = c('leiden_clus', 'cell_types', 'nr_genes'))

annotated_grid_metadata

mini_visium = createSpatialNetwork(gobject = mini_visium, minimum_k = 2, maximum_distance_delaunay = 400)
mini_visium = createSpatialNetwork(gobject = mini_visium, minimum_k = 2, method = 'kNN', k = 10)
showNetworks(mini_visium)

p1 <- spatPlot(gobject = mini_visium, show_network = T,
               network_color = 'blue', spatial_network_name = 'Delaunay_network',
               point_size = 2.5, cell_color = 'leiden_clus') + ggtitle("Delaunay_network")

p2 <- spatPlot(gobject = mini_visium, show_network = T,
               network_color = 'blue', spatial_network_name = 'kNN_network',
               point_size = 2.5, cell_color = 'leiden_clus') + ggtitle("kNN_network")

p1 + p2

km_spatialgenes = binSpect(mini_visium)
km_spatialgenes

spatGenePlot(mini_visium, expression_values = 'scaled',
             genes = km_spatialgenes[1:4]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

rank_spatialgenes = binSpect(mini_visium, bin_method = 'rank')
spatGenePlot(mini_visium, expression_values = 'scaled',
             genes = rank_spatialgenes[1:4]$genes,
             point_shape = 'border', point_border_stroke = 0.1,
             show_network = F, network_color = 'lightgrey', point_size = 2.5,
             cow_n_col = 2)

ext_spatial_genes = km_spatialgenes[1:100]$genes
ext_spatial_genes

# 检测空间相关的基因
spat_cor_netw_DT = detectSpatialCorGenes(mini_visium,
                                         method = 'network', spatial_network_name = 'Delaunay_network',
                                         subset_genes = ext_spatial_genes)

# 2. cluster correlation scores
spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)
heatmSpatialCorGenes(mini_visium, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus')

# 可视化和过滤空间相关的基因
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 6, show_top_genes = 1)
top_netw_spat_cluster

cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
cluster_genes
mini_visium = createMetagenes(mini_visium, gene_clusters = cluster_genes, name = 'cluster_metagene')
spatCellPlot(mini_visium,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 3,show_image=T,show_legend=F)

hmrf_folder = paste0(temp_dir,'11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

?doHMRF
# perform hmrf
my_spatial_genes = km_spatialgenes[1:100]$genes
HMRF_spatial_genes = doHMRF(gobject = mini_visium,
                            expression_values = 'scaled',
                            spatial_genes = my_spatial_genes,
                            spatial_network_name = 'Delaunay_network',
                            k = 8,
                            betas = c(28,2,2),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes_brain/SG_top100_k8_scaled'))

# check and select hmrf
myviewHMRFresults2D <- edit(viewHMRFresults2D)
environment(myviewHMRFresults2D) <- environment(viewHMRFresults2D)
debug(myviewHMRFresults2D)
?viewHMRFresults2D
viewHMRFresults2D(gobject = mini_visium,
                  HMRFoutput = HMRF_spatial_genes,
                  k = 8, betas_to_view = 2,
                  point_size = 2)

for(i in seq(28, 30, by = 2)) {
  myviewHMRFresults2D(gobject = mini_visium,
                      HMRFoutput = HMRF_spatial_genes,
                      k = 8, betas_to_view = i,
                      point_size = 2)
}

myaddHMRF <- edit(addHMRF)
environment(myaddHMRF) <- environment(addHMRF)
mini_visium = myaddHMRF(gobject = mini_visium,
                        HMRFoutput = HMRF_spatial_genes,
                        k = 8, betas_to_add = c(28),
                        hmrf_name = 'HMRF')

mini_visium@cell_metadata
giotto_colors = getDistinctColors(8)
giotto_colors
names(giotto_colors) = 1:8
head(mini_visium@cell_metadata)
spatPlot(gobject = mini_visium, cell_color = 'HMRF_k8_b.28',
         point_size = 3, coord_fix_ratio = 1, cell_color_code = giotto_colors)

