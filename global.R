if(! "librarian" %in% installed.packages()){
  install.packages("librarian", dependencies = TRUE)
}
if(! "CIPR" %in% installed.packages()){
  devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = F)
}



require("librarian")
librarian::shelf(shinydashboard, NeuCA, DT, shiny, shinyjs, shinycssloaders, shinyalert, SeuratObject, Seurat, plotly, igraph,
                 hdf5r, rgl, RColorBrewer, dplyr, visNetwork, heatmaply, gprofiler2, ggplot2, ggpubr, CIPR, dittoSeq, slingshot,
                 saeyslab/nichenetr,tidyverse,theislab/destiny,carmonalab/UCell, colorspace, missMDA, dismo, scDblFinder,phateR,
                 decoupleR, tibble, tidyr, shinyBS, glmGamPoi, HGNChelper, openxlsx, data.tree, ggraph, shinyChakraSlider)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Kidney" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
#ATAC libraries
# ArchR
# pheatmap
# SCENIC
# SCopeLoomR
# AUCell
# GSEABase
# RcisTarget
# stringr
# readr
# parallel
# chromVAR
# chromVARmotifs
# reticulate
# JASPAR2020
# JASPAR2018
# JASPAR2016 #BiocManager::install("JASPAR2020", BiocManager::install("JASPAR2018", BiocManager::install("JASPAR2016"

#Global variables
YEAR <- substr(Sys.Date(), 1, 4)
user_dir <- "" #user's folder in temp
user_dir_pyscenic <- "" #user's folder for pyscenic results

#tab Upload
seurat_object <- NULL 
init_seurat_object <- NULL
#my_metadata <- NULL
minimum_cells <- 3
minimum_features <- 200
organism <- "mouse" #or human

#tab Quality Control
qc_minFeatures <- 500
qc_maxFeatures <- 6000
qc_maxMtPercent <- 10

#tab Normalization
normalize_normMethod <- "LogNormalize"
normalize_normScaleFactor <- 10000
normalize_hvgMethod <- "vst"
normalize_hvgNGenes <- 2000
normalize_scaleRegressOut <- NULL

#tab Clustering
snn_dims <- 15
snn_k <- 20
cluster_res <- 0.6
cluster_dims <- 15
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

#tab DEA
markers_test <- "wilcox"
markers_minPct <- "0.1"
markers_minLogfc <- "0.25"
markers_minPval <- "0.01"
markers_logFCBase <- "avg_logFC"

#tabs Umap/tsne, DEA, Cell cycle, Trajectory
reductions_choices <- c("-")

#Cell cycle
mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
#export tables RNA

export_metadata_RNA <- ""
export_loadingScoresTable_RNA <- ""
export_clustertable_RNA <- ""
export_markerGenes_RNA <- ""
export_enrichedTerms_RNA <- ""
export_annotation_RNA <- ""
export_ligandReceptor_full_RNA <- ""
export_ligandReceptor_short_RNA <- ""
export_scenicAUC_full_RNA <- ""
export_decoupleRZscores_full_RNA <- ""

#ATAC variables
ArrowFiles <- NULL
proj_default <- NULL
#export tables
export_metadata_ATAC <- ""
export_clustertable_ATAC <- ""
export_markerGenes_ATAC <- ""
export_markerPeaks_ATAC <- ""
export_motifs_ATAC <- ""
export_positiveRegulators_ATAC <- ""
export_peakToGenelinks_ATAC <- ""
export_PeakMotifTable_ATAC <- ""

userMode <- FALSE

if(userMode == F)
{
  #to improve speed
}

#opening new window
js.enrich <- "
  shinyjs.Enrich = function(url) {
    window.open(url[0]);
  }
"

# Functions ####

# This is a void function that hides all shiny css loaders.
hideAllLoaders <- function(){
  shinyjs::hide("hvgScatter_loader")
  shinyjs::hide("nFeatureViolin_loader")
  shinyjs::hide("totalCountsViolin_loader")
  shinyjs::hide("mitoViolin_loader")
  shinyjs::hide("genesCounts_loader")
  shinyjs::hide("mtCounts_loader")
  shinyjs::hide("filteredNFeatureViolin_loader")
  shinyjs::hide("filteredTotalCountsViolin_loader")
  shinyjs::hide("filteredMitoViolin_loader")
  shinyjs::hide("filteredGenesCounts_loader")
  shinyjs::hide("filteredMtCounts_loader")
  shinyjs::hide("TSS_plot_loader")
  shinyjs::hide("nFrag_plot_loader")
  shinyjs::hide("TSS_nFrag_plot_loader")
  shinyjs::hide("elbowPlotPCA_loader")
  shinyjs::hide("PCAscatter_loader")
  shinyjs::hide("PCAloadings_loader")
  shinyjs::hide("PCAheatmap_loader")
  shinyjs::hide("clusterBarplot_loader")
  shinyjs::hide("clusterBarplotATAC_loader")
  shinyjs::hide("umapPlot_loader")
  shinyjs::hide("umapPlotATAC_loader")
  shinyjs::hide("findMarkersHeatmap_loader")
  shinyjs::hide("findMarkersDotplot_loader")
  shinyjs::hide("findMarkersFeaturePlot_loader")
  shinyjs::hide("findMarkersFPfeature1_loader")
  shinyjs::hide("findMarkersFPfeature2_loader")
  shinyjs::hide("findMarkersFPfeature1_2_loader")
  shinyjs::hide("findMarkersFPcolorbox_loader")
  shinyjs::hide("findMarkersViolinPlot_loader")
  shinyjs::hide("findMarkersVolcanoPlot_loader")
  shinyjs::hide("findMarkersFeaturePlotATAC_loader")
  shinyjs::hide("snnSNN_loader")
  shinyjs::hide("findMarkersGenesHeatmapATAC_loader")
  shinyjs::hide("findMarkersGenesATACTable_loader")
  shinyjs::hide("findMarkersPeaksATACTable_loader")
  shinyjs::hide("findMarkersPeaksHeatmapATAC_loader")
  shinyjs::hide("doubletATAC_loader3")
  shinyjs::hide("doubletATAC_loader4")
  shinyjs::hide("cellCyclePCA_loader")
  shinyjs::hide("cellCycleBarplot_loader")
  shinyjs::hide("gProfilerManhattan_loader")
  shinyjs::hide("findMotifsHeatmapATAC_loader")
  shinyjs::hide("findMotifsATACTable_loader")
  shinyjs::hide("annotateClustersCIPRDotplot_loader")
  shinyjs::hide("annotateClustersUMAP_loader")
  shinyjs::hide("ligandReceptorFullHeatmap_loader")
  shinyjs::hide("ligandReceptorCuratedHeatmap_loader")
  shinyjs::hide("trajectoryPlot_loader")
  shinyjs::hide("trajectoryPseudotimePlot_loader")
  shinyjs::hide("trajectoryPseudotimePlotATAC_loader")
  shinyjs::hide("grnHeatmapRNA_loader")
  shinyjs::hide("grnHeatmapATAC_loader")
  shinyjs::hide("grnATACTable_loader")
  shinyjs::hide("grnATACTable2_loader")
  shinyjs::hide("grnATACTable3_loader")
  shinyjs::hide("visualizeTracksOutput_loader")
}

