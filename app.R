library(bsplus)
library(shiny)
library(Seurat)
library(Signac)
library(shinydashboard)
library(shinyjs)
library(plotly)
source("global.R", local=TRUE)
js.enrich <- "
  shinyjs.Enrich = function(url) {
    window.open(url[0]);
  }"
shinyApp(
  
  ui = dashboardPage(
    title="KFO single-cell analysis",
    skin = "purple",
    #Header
    dashboardHeader(
      title = "Single-cell analysis"
    ),
    #Sidebar
    dashboardSidebar(sidebarMenu(id = "sidebarMenu",
                                 menuItem(text = "Main page", tabName = "home", icon = icon("home")),
                                 tags$hr(),
                                 menuItem(text = "1. Data upload", tabName = "input"),
                                 menuItem(text = "2. Droplet processing", tabName = "dp"),
                                 menuItem(text = "3. Quality control", tabName = "qc"),
                                 menuItem(text= "4. Data normalization", tabName = "scnormalize"),
                                 menuItem(text = "5. PCA", tabName = "pca"),
                                 menuItem(text = "6. Clustering", tabName = "cluster"),
                                 menuItem(text = "7. UMAP", tabName = "umap"),
                                 menuItem(text = "8. Cell cycle analysis", tabName = "cellcycle"),
                                 menuItem(text = "9. Dataset inspection", tabName = "shinycell"),
                                 menuItem(text = "10. Marker gene detection", tabName = "markers"),
                                 menuItem(text = "11. Cell-type annotation", tabName = "annotateClusters"),
                                 tags$hr(),
                                 menuItem(text= "Utilities",tabName = "utilities")
                                # menuItem(text = "Help", tabName = "help", icon = icon("question"))
    )),
    #Body
    dashboardBody(tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "main.css")),
                  tags$head(tags$link(rel = "stylesheet", type = "text/css", href="loading-bar.css")), # loading bar CSS
                  tags$head(tags$script(src = "rshiny_handlers.js")), # R to JS
                  tags$head(tags$script(src = "loading-bar.js")), # loading bar JS
                  useShinyjs(),
                  extendShinyjs(text = js.enrich, functions = c("Enrich")),
                  tabItems(
      tabItem(tabName = "home",
              div(id = "home_div", class = "div_container",
                  h1(class = "container_title", "KFO 329 single-cell Shiny app"),
                  HTML("<p class=container_text> This application is made to perform the standard pipeline for analyzing scRNA-seq datasets in the KFO. 
                       You can upload your data or try the example sample."))),
      tabItem(tabName = "input",
              fluidRow(
                box(width = 3, status = "info", solidHeader = TRUE,
                           title = "Data upload",
                           tabsetPanel(type="tabs",
                                       tabPanel("10x input files (MEX)",
                                                tags$h3("File upload"),
                                                tags$hr(),
                                                textInput(inputId = "upload10xRNAprojectID", label = "Project name : ", value = "KFO 329"),
                                                fileInput(inputId = "barcodes", label = "1. Choose barcodes.tsv.gz file", accept = ".gz"),
                                                fileInput(inputId = "genes", label = "2. Choose features.tsv.gz file", accept = ".gz"),
                                                fileInput(inputId = "matrix", label = "3. Choose matrix.mtx.gz file", accept = ".gz"),
                                                sliderInput(inputId = "upload10xRNAminCells", label = "Include genes detected in at least this many cells :", min = 0, max = 20, value = 3, step = 1),
                                                sliderInput(inputId = "upload10xRNAminFeatures", label = "Include cells where at least this many genes are detected :", min = 0, max = 1000, value = 200, step = 1),
                                                radioButtons("upload10xRNARadioSpecies", label = h3("Select organism : "),
                                                             choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                            "Homo sapiens (Human)" = "human"
                                                             ), 
                                                             selected = "mouse"),
                                                actionButton(inputId = "upload10xRNAConfirm", label = "Submit", class="btn-run", icon = icon("check-circle")),
                                                tags$h3("Export working object as .RDS file"),
                                                tags$hr(),
                                                downloadButton(outputId = "utilitiesConfirmExport2", label = "Export .RDS")
                                       ),
                                       tabPanel("10x input files (HDF5)",
                                                tags$h3("File upload"),
                                                tags$hr(),
                                                textInput(inputId = "uploadh5projectID", label = "Project name : ", value = "KFO 329"),
                                                fileInput(inputId = "uploadh5", label = "1. HDF5 file", accept = ".h5"),
                                                sliderInput(inputId = "uploadh5minCells", label = "Include genes detected in at least this many cells :", min = 1, max = 20, value = 3, step = 1),
                                                sliderInput(inputId = "uploadh5minFeatures", label = "Include cells where at least this many genes are detected :", min = 1, max = 1000, value = 200, step = 1),
                                                radioButtons("uploadh5RadioSpecies", label = h3("Select organism : "),
                                                             choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                            "Homo sapiens (Human)" = "human"
                                                             ), 
                                                             selected = "mouse"),
                                                actionButton(inputId = "uploadh5Confirm", label = "Submit", class="btn-run", icon = icon("check-circle")),
                                                tags$h3("Export working object as .RDS file"),
                                                tags$hr(),
                                                downloadButton(outputId = "utilitiesConfirmExport1", label = "Export .RDS")
                                       ),
                                       tabPanel("Seurat object (.rds)",
                                                fileInput(inputId = "uploadRdsFile", label = "Choose a Seurat object saved in .rds format", accept = ".rds"),
                                                radioButtons("uploadRdsRadioSpecies", label = h3("Select organism : "),
                                                             choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                            "Homo sapiens (Human)" = "human"
                                                             ), 
                                                             selected = "mouse"),
                                                actionButton(inputId = "uploadSeuratRdsConfirm", label = "Load Seurat object", class="btn-run", icon = icon("check-circle")),
                                                tags$br(),
                                                tags$br(),
                                                tags$hr(),
                                                tags$br(),
                                                selectInput("utilitiesActiveAssay", "(Optional) Change active assay:",
                                                            c("Assay" = "RNA")),
                                                actionButton(inputId = "utilitiesConfirmChangeAssay", label = "Change assay"),
                                                tags$br(),
                                                tags$h3("Export working object as .RDS file"),
                                                tags$hr(),
                                                downloadButton(outputId = "utilitiesConfirmExport", label = "Export .RDS")
                                                )
                                       )
                           ),
                box( width = 8, solidHeader = TRUE, status = "info",
                     title = "Metadata table",
                     #tabsetPanel(type = "tabs", id = "uploadTabPanel",
                      #           tabPanel("scRNA-seq",
                              dataTableOutput("metadataTable"),
                              downloadButton(outputId = "uploadMetadataExportRNA", label = "Save table")
                                # )
                    # )
                     )
                       )
              
              ),
      
      #utilities tab
      tabItem(tabName = "utilities",
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Edit/download Seurat object",
                  tags$h3("Active clusters"),
                  tags$hr(),
                  selectInput("utilitiesActiveClusters", "Set active clustering column:",
                              c("orig.ident" = "orig.ident")),
                  actionButton(inputId = "utilitiesConfirmChangeCluster", label = "Change clustering column", class="btn-run", icon = icon("check-circle")),
                  tags$h3("Rename cluster"),
                  tags$hr(),
                  selectInput(inputId = "utilitiesRenameOldName", label = "Cluster to be renamed (old name):", choices = "-", multiple = F),
                  textInput(inputId = "utilitiesRenameNewName", label = "New cluster name:", value = "New_name"),
                  actionButton(inputId = "utilitiesConfirmRename", label = "Rename", class="btn-run", icon = icon("check-circle")),
                  tags$h3("Delete cluster"),
                  tags$hr(),
                  selectInput(inputId = "utilitiesDeleteCluster", label = "Cluster to be deleted:", choices = "-", multiple = F),
                  actionButton(inputId = "utilitiesConfirmDelete", label = "Delete", class="btn-run", icon = icon("check-circle")),
                  tags$h3("Command history"),
                  tags$hr(),
                  actionButton(inputId = "utilitiesPlotCommands", label = "View command history", class="btn-run", icon = icon("check-circle")),
                  verbatimTextOutput(outputId = "history")
                )
              )
      ),
      #QC tab
      tabItem(tabName = "qc",
              tabsetPanel(type = "tabs", id = "qcTabPanel",
                          fluidRow(
                            box(
                              width = 3, status = "info", solidHeader = TRUE,
                              title = "Quality control",
                              tags$h3("1. Display quality control plots before filtering"),
                              actionButton(inputId = "qcDisplay", label = "Display plots", class="btn-run", icon = icon("check-circle")),
                              tags$hr(),
                              tags$h3("2. Filter out low quality cells"),
                              tags$hr(),
                              sliderInput(inputId = "minUniqueGenes", label = "Minimum genes detected", min = 200, max = 2000, value = 500, step = 1)%>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Filter out cells that have unique gene counts less than:", placement = "left"
                                    )
                                ),
                              sliderInput(inputId = "maxUniqueGenes", label = "Maximum genes detected", min = 2001, max = 7000, value = 4500, step = 1)%>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Filter out cells that have unique gene counts over than:", placement = "left"
                                    )
                                ),
                              
                              sliderInput(inputId = "maxMtReads", label = "Mitochondrial percentage", min = 0.1, max = 10, value = 1, step = 0.1)%>%
                                shinyInput_label_embed(
                                  shiny_iconlink() %>%
                                    bs_embed_popover(
                                      title = "Filter out cells that their percentage of genes mapped to mitochondrial genome exceeds:", placement = "left"
                                    )
                                ),
                              selectInput("qcColorBy", "Color by:",
                                          c("orig.ident" = "orig.ident")),
                              actionButton(inputId = "qcConfirm", label = "Perform filtering", class="btn-run", icon = icon("check-circle"))
                            ),
                            box(
                              width = 9, status = "info", solidHeader = TRUE,
                              title = "Quality control plots",
                              
                              tabsetPanel(type="tabs", id = "qc_tabs_rna",
                                          tabPanel("Pre-filtering plots",
                                                   column(
                                                     div(id="nFeatureViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "nFeatureViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="totalCountsViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "totalCountsViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="mitoViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "mitoViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="genesCounts_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "genesCounts", height= "100%")
                                                         )
                                                     ), width = 6),
                                                   column(
                                                     div(id="mtCounts_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "mtCounts", height= "100%")
                                                         )
                                                     ), width = 6),
                                                   column(verbatimTextOutput(outputId = "cellStats"), width = 4)
                                          ),
                                          tabPanel("Post-filtering plots",
                                                   column(
                                                     div(id="filteredNFeatureViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "filteredNFeatureViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="filteredTotalCountsViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "filteredTotalCountsViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="filteredMitoViolin_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "filteredMitoViolin", height = "100%")
                                                         )
                                                     ), width = 4),
                                                   column(
                                                     div(id="filteredGenesCounts_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "filteredGenesCounts", height= "100%")
                                                         )
                                                     ), width = 6),
                                                   column(
                                                     div(id="filteredMtCounts_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId = "filteredMtCounts", height= "100%")
                                                         )
                                                     ), width = 6),
                                                   column(verbatimTextOutput(outputId = "filteredCellStats"), width = 4)
                                          )
                              
                              )
                            )
                          )
                          )
      ),
      
#Droplet

tabItem(tabName = "dp", 
        tags$br(),
        tabsetPanel(type = "tabs", id = "doubletDetectionTabPanel",
                    tabPanel("Doublet detection",
                             fluidRow(
                               box(
                                 width = 4, status = "info", solidHeader = TRUE,
                                 title = "Doublet detection parameters",
                                 tags$h3("1. Options for doublet detection"),
                                 tags$hr(),
                                 textInput(inputId = "doubletsPN", label = "Expected doublet rate:", value=NULL),
                                 textInput(inputId = "doubletsart", label = "Number of genes to use:", value=NULL),
                               selectInput(inputId = "samplecolumn", label="Sample column:", choices =  c("orig.ident" = "orig.ident")),
                                 actionButton(inputId = "doubletsConfirm", label = "Perform doublet detection", class="btn-run", icon = icon("check-circle")),
                                 tags$h3("2. Doublet removal (optional)"),
                                 actionButton(inputId = "doubletsRemove", label = "Remove doublets", class="btn-run", icon = icon("check-circle"))
                               ),
                               box(
                                 width = 8, status = "info", solidHeader = TRUE, title = "Doublet detection output",
                                 verbatimTextOutput(outputId = "doubletsInfo")
                               )
                             )
                    )
        )
),
#Normalization
      tabItem(tabName = "scnormalize",
              tags$br(),
              fluidRow(
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "Apply sctransform normalization",
                  tags$hr(),
                  selectInput("normalizeRegressColumns", "Select variables to regress out", list(), selected = NULL, multiple = TRUE, selectize = TRUE, width = NULL, size = NULL),
                  actionButton(inputId = "normalizeConfirm", label = "Run sctransform", class="btn-run", icon = icon("check-circle"))
                ),
                box(
                  width = 8, status = "info", solidHeader = TRUE,
                  title = "Highly variable genes",
                  div(id="hvgScatter_loader",
                      shinycssloaders::withSpinner(
                        plotlyOutput(outputId = "hvgScatter", height = "800px")
                      )
                  ),
                  column(verbatimTextOutput(outputId = "hvgTop10Stats"), width = 8)
                )
              )
      ),
#PCA tab
tabItem(tabName = "pca",
        tags$br(),
        tabsetPanel(type = "tabs", id = "pcaTabPanel",
                    fluidRow(
                      box(
                        width = 12, status = "info", solidHeader = TRUE,
                        title = "PCA results", height = "1200px",
                        tabsetPanel(type = "tabs",
                                    tabPanel("PCA run",
                                             column(actionButton(inputId = "PCrunPCA", label = "Run PCA", class="btn-run", icon = icon("check-circle")), width = 12),
                                             div(
                                               column(
                                                 div(id="elbowPlotPCA_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "elbowPlotPCA", height = "750px")
                                                     )
                                                 ), width = 6),
                                               column(
                                                 div(id="PCAscatter_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "PCAscatter", height = "750px")
                                                     )
                                                 ), width = 6)
                                             )
                                    ),
                                    tabPanel("PCA exploration",
                                             selectInput("PCin", "Select a principal component :", choices=1:100, selected = 1, multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                                             column(actionButton(inputId = "PCconfirm", label = "Explore the principal component", class="btn-run", icon = icon("check-circle")),
                                                    width = 12),
                                             div(
                                               column(
                                                 tags$h3("PCA loading scores (top-30 genes for this PC)"),
                                                 div(id="PCAloadings_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "PCAloadings", height = "700px")
                                                     )
                                                 ), width = 6),
                                               column(
                                                 tags$h3("Heatmap of scaled expression (top-30 genes for this PC)"),
                                                 div(id="PCAheatmap_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "PCAheatmap", height = "700px")
                                                     )
                                                 ), width = 6)
                                             ),
                                             downloadButton(outputId = "pcaRNAExport", label = "Save table")
                                    )
                        )
                        )
                    )
        )
        
        ),

#Clustering tab
tabItem(tabName = "cluster", 
        tags$br(),
        tabsetPanel(type = "tabs", id = "clusteringTabPanel",
                    fluidRow(
                      box(
                        width = 4, status = "info", solidHeader = TRUE,
                        title = "Clustering options",
                        tags$h3("1. Construction of the shared nearest neighbour (SNN) graph"),
                        tags$hr(),
                        sliderInput(inputId = "snnK", label = "Number of neighbours for each cell [k]:", min = 1, max = 200, value = 20, step = 1),
                        sliderInput(inputId = "snnPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                        tags$h3("2. Communities' detection (Louvain algorithm)"),
                        tags$hr(),
                        sliderInput(inputId = "clusterRes", label = "Clustering resolution :", min = 0.1, max = 5, value = 0.5, step = 0.1),
                        actionButton(inputId = "snnConfirm", label = "Run clustering", class="btn-run", icon = icon("check-circle"))
                      ),
                      box(
                        width = 8, status = "info", solidHeader = TRUE, title = "Clustering output",
                        tabsetPanel(type = "tabs",
                                    tabPanel("Clustering results",
                                             tabsetPanel(type = "tabs",
                                                         tabPanel("Cluster table",
                                                                  dataTableOutput(outputId="clusterTable"),
                                                                  downloadButton(outputId = "clusterTableRNAExport", label = "Save table")
                                                         ),
                                                         tabPanel("Cluster barplot",
                                                                  selectInput("clusterGroupBy", "Grouping variable:",
                                                                              c("orig.ident" = "orig.ident")),
                                                                  actionButton(inputId = "clusterBarplotConfirm", label = "Display barchart", class="btn-run", icon = icon("check-circle")),
                                                                  div(id="clusterBarplot_loader",
                                                                      shinycssloaders::withSpinner(
                                                                        plotlyOutput(outputId = "clusterBarplot", height = "700px")
                                                                      )
                                                                  )
                                                         )
                        )
                      ),
                      tabPanel("Shared Nearest Neighbour (SNN) Graph", 
                               actionButton(inputId = "snnDisplayConfirm", label = "Display SNN graph", class="btn-run", icon = icon("check-circle")),
                               div(id="snnSNN_loader",
                                   shinycssloaders::withSpinner(
                                     visNetworkOutput(outputId="snnSNN", height = "1300px")
                                   )
                               )
                      )
                    )
        )
        
)
)
),
#Marker analysis
tabItem(tabName = "markers", 
        tags$br(),
        tags$br(),
        tabsetPanel(type = "tabs", id = "markers",
                    tabPanel("scRNA-seq",
                             fluidRow(
                               box(width = 3, status = "info", solidHeader = TRUE,
                                   title = "Differential Expression Analysis options", 
                                   selectInput("findMarkersTest", "Test used:",
                                               c("Wilcoxon rank sum test" = "wilcox",
                                                 "Likelihood-ratio test for single cell feature expression" = "bimod",
                                                 "Standard AUC classifier" = "roc",
                                                 "Student's t-test" = "t",
                                                 "MAST" = "MAST",
                                                 "DESeq2" = "DESeq2"
                                               )),
                                   radioButtons("findMarkersLogBase", label = "Base used for average logFC calculation: ",
                                                choices = list("log(e)" = "avg_logFC", 
                                                               "log(2)" = "avg_log2FC"
                                                ), 
                                                selected = "avg_log2FC"),
                                   
                                   sliderInput(inputId = "findMarkersMinPct", label = "Minimum % of expression", min = 0, max = 1, value = 0.1, step = 0.01)%>%
                                     shinyInput_label_embed(
                                       shiny_iconlink() %>%
                                         bs_embed_popover(
                                           title = "Only test genes that are detected in a minimum fraction of cells in either of the two populations:", placement = "bottom"
                                         )
                                     ),
                                   
                                   sliderInput(inputId = "findMarkersLogFC", label = "Avg Log FC threshold", min = 0, max = 3, value = 0.25, step = 0.05)%>%
                                     shinyInput_label_embed(
                                       shiny_iconlink() %>%
                                         bs_embed_popover(
                                           title = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells:", placement = "bottom"
                                         )
                                     ),
                                   
                                   sliderInput(inputId = "findMarkersPval", label = "P-value threshold", min = 0, max = 1, value = 0.05, step = 0.01)%>%
                                     shinyInput_label_embed(
                                       shiny_iconlink() %>%
                                         bs_embed_popover(
                                           title = "Only return markers that have a p-value < slected threshold, or a power > selected threshold (if the test is ROC) :", placement = "bottom"
                                         )
                                     ),
                                   actionButton(inputId = "findMarkersConfirm", label = "Find marker genes", class="btn-run", icon = icon("check-circle"))
                               ),
                               
                               box(
                                 width = 9, status = "info", solidHeader = TRUE, title = "DEA results",
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Marker genes", 
                                                      dataTableOutput(outputId="findMarkersTable"),
                                                      downloadButton(outputId = "findMarkersRNAExport", label = "Save table")),
                                             tabPanel("Heatmap", 
                                                      actionButton(inputId = "findMarkersTop10HeatmapConfirm", label = "Display top-10 marker genes heatmap", class="btn-run", icon = icon("check-circle")),
                                                      div(id="findMarkersHeatmap_loader",
                                                          shinycssloaders::withSpinner(
                                                            plotlyOutput(outputId = "findMarkersHeatmap", height = "1300px")
                                                          )
                                                      )
                                             ),
                                             
                                             tabPanel("Dotplot", 
                                                      actionButton(inputId = "findMarkersTop10DotplotConfirm", label = "Display top-10 marker genes dotplot", class="btn-run", icon = icon("check-circle")),
                                                      div(id="findMarkersDotplot_loader",
                                                          shinycssloaders::withSpinner(
                                                            plotlyOutput(outputId = "findMarkersDotplot", height = "1300px")
                                                          )
                                                      )
                                             ),
                                             tabPanel("VolcanoPlot", fluidRow(
                                               box(width = 3, status = "info", solidHeader = TRUE, title = "Cluster selection",
                                                   selectInput("findMarkersClusterSelect", "Cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                                                   actionButton(inputId = "findMarkersVolcanoConfirm", "Display volcano plot", class="btn-run", icon = icon("check-circle"))
                                               ),
                                               
                                               box(width = 9, status = "info", solidHeader = TRUE, title = "Volcano plot",
                                                   div(id="findMarkersVolcanoPlot_loader",
                                                       shinycssloaders::withSpinner(
                                                         plotlyOutput(outputId = "findMarkersVolcanoPlot", height = "800px")
                                                       )
                                                   )
                                               )
                                             )
                                             )
                                 )
                               )
                             )
                    )
        )	
),

#UMAP tab
tabItem(tabName = "umap", 
        tags$br(),
       tabsetPanel(type = "tabs", id = "umapTabPanel",
                    tabPanel("scRNA-seq",
                             fluidRow(
                               box(width = 3, status = "info", solidHeader = TRUE,
                                   title = "Cells visualization options",
                                   sliderInput(inputId = "umapSeed", label = "Set seed :", min = 1, max = 500, value = 42, step = 1),
                                   sliderInput(inputId = "umapPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                   actionButton(inputId = "umapRunUmap", label = "Run UMAP", class="btn-run", icon = icon("check-circle")),
                                   actionButton(inputId = "umapRunTsne", label = "Run tSNE", class="btn-run", icon = icon("check-circle")),
                                   actionButton(inputId = "umapRunDFM", label = "Run Diffusion Map", class="btn-run", icon = icon("check-circle")),
                                   actionButton(inputId = "umapRunPhate", label = "Run PHATE", class="btn-run", icon = icon("check-circle")),
                                   tags$h3("Display settings"),
                                   tags$hr(),
                                   selectInput("umapType", "Plot type:",
                                               c("-" = "-")
                                   ),
                                   selectInput("umapDimensions", "Dimensions:",
                                               c("2D" = "2",
                                                 "3D" = "3")),
                                   selectInput("umapColorBy", "Color by:",
                                               c("Cluster" = "seurat_clusters")),
                                   
                                   sliderInput("umapDotSize", "Size:", min = 1, max = 20, value = 5, step = 0.5),
                                   sliderInput("umapDotOpacity", "Opacity:", min = 0, max = 1, value = 1, step = 0.1),
                                   sliderInput("umapDotBorder", "Border width:", min = 0, max = 10, value = 0.5, step = 0.1),
                                   actionButton(inputId = "umapConfirm", label = "Update plot", class="btn-run", icon = icon("check-circle"))
                               ),
                               
                               box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                                   div(id="umapPlot_loader",
                                       shinycssloaders::withSpinner(
                                         plotlyOutput(outputId = "umapPlot", height = "1100px")
                                       )
                                   )
                               )
                             )
                    )
        )
),
#Cell cycle phase analysis
tabItem(tabName = "cellcycle",
        tags$br(),
        fluidRow(
          box(
            width = 12, status = "info", solidHeader = T,
            title = "Cell cycle analysis",
            tabsetPanel(type = "tabs",
                        tabPanel("Dimensionality reduction plot", 
                                 tags$br(),
                                 tags$div("After succesfully running cell cycle 
                                       phase analysis the columns S.Score, G2M.Score and CC.Difference are stored in the metadata
                                       table. You can regress out the cell cycle effect by returning to the tab data normalization
                                          and repeating step 4.",id="cellCycleMessage"),
                                 tags$br(),
                                 tags$hr(),
                                 selectInput("cellCycleReduction", "Plot type:",
                                             c("-" = "-")
                                 ),
                                 actionButton(inputId = "cellCycleRun", label = "Cell cycle analysis", class="btn-run", icon = icon("check-circle")),
                                 div(id="cellCyclePCA_loader",
                                     shinycssloaders::withSpinner(
                                       plotlyOutput(outputId = "cellCyclePCA", height = "700px")
                                     )
                                 )
                        ),
                        tabPanel("Barplot",
                                 div(id="cellCycleBarplot_loader",
                                     shinycssloaders::withSpinner(
                                       plotlyOutput(outputId = "cellCycleBarplot", height = "3100px")
                                     )
                                 )
                        )
            )
          )
        )
),
#Clusters' annotation
tabItem(tabName = "annotateClusters",
        tags$br(),
        tabsetPanel(type="tabs", id = "annotateClustersTabPanel",
                    tabPanel("ScType",
                             fluidRow(
                             box(
                               width = 3, status = "info", solidHeader = TRUE,
                               actionButton(inputId = "annotateClustersConfirmtype", label = "Run cluster annotation analysis", class="btn-run", icon = icon("check-circle"))
                               
                             ),
                             box(
                               width = 9, status = "info", solidHeader = TRUE, title = "scType plots",
                               tabsetPanel(type = "tabs",
                                           tabPanel("Bubble Plot", 
                                                    plotlyOutput(outputId="annotateClustersBubbleplot", height = "1100px")
                                           ),
                                           tabPanel("Tissue type prediction", 
                                                    div(id="annotateClustersCIPRDotplot_loader",
                                                        shinycssloaders::withSpinner(
                                                          plotlyOutput(outputId="annotateClustersTissueplot", height = "1100px")
                                                        )
                                                    )
                                           )
                               ))
                             )
                             ),
                    tabPanel("CIPR",
                             fluidRow(
                               box(
                                 width = 3, status = "info", solidHeader = TRUE,
                                 title = "Annotation parameters",
                                 radioButtons("annotateClustersReference", label = "Reference dataset : ",
                                              choices = list("ImmGen (mouse)" = "immgen", 
                                                             "Presorted RNAseq (mouse)" = "mmrnaseq",
                                                             "Blueprint-Encode (human)" = "blueprint",
                                                             "Primary Cell Atlas (human)" = "hpca",
                                                             "DICE (human)" = "dice",
                                                             "Hematopoietic diff (human)" = "hema",
                                                             "Presorted RNA seq (human)" = "hsrnaseq"
                                              ),
                                              selected = "mmrnaseq"
                                 ),
                                 tags$hr(),
                                 sliderInput("annotateClustersSlider", "Keep top Nth % of variable genes in reference :", min = 0, max = 100, value = 100, step = 1),
                                 tags$hr(),
                                 radioButtons("annotateClustersMethod", label = "Select method for comparisons : ",
                                              choices = list("logFC dot product" = "logfc_dot_product", 
                                                             "logFC Spearman" = "logfc_spearman",
                                                             "logFC Pearson" = "logfc_pearson",
                                                             "Spearman (all genes)" = "all_genes_spearman",
                                                             "Pearson (all genes)" = "all_genes_pearson"
                                              ),
                                              selected = "all_genes_pearson"
                                 ),
                                 tags$hr(),
                                 actionButton(inputId = "annotateClustersConfirm", label = "Run cluster annotation analysis", class="btn-run", icon = icon("check-circle")),
                               ),
                               box(
                                 width = 9, status = "info", solidHeader = TRUE, title = "Cell type annotation",
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Top-5 hits table", 
                                                      dataTableOutput(outputId="annotateClustersCIPRTable"),
                                                      downloadButton(outputId = "annotationRNAExport", label = "Save table")
                                             ),
                                             tabPanel("Top-5 hits dotplot", 
                                                      div(id="annotateClustersCIPRDotplot_loader",
                                                          shinycssloaders::withSpinner(
                                                            plotlyOutput(outputId="annotateClustersCIPRDotplot", height = "1100px")
                                                          )
                                                      )
                                             )
                                 ))
                                                )
        ) )
),
#Feature inspection
tabItem(tabName = "shinycell",
        tags$br(),
        tabsetPanel(type = "tabs", id = "featuresTabPanel",
                    tabPanel("scRNA-seq",
                             fluidRow(
                               tabsetPanel(type = "tabs",
                                           tabPanel("Feature plot", fluidRow(
                                             box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                 selectizeInput(inputId = 'findMarkersGeneSelect',
                                                                label = 'Type the name of a gene, gene signature or numeric metadata column: 
                                                                      (e.g. "Ccl2" or "Signature1_Ucell", or "nCount_RNA")',
                                                                choices = NULL,
                                                                selected = NULL,
                                                                multiple = FALSE),
                                                 selectInput("findMarkersReductionType", "Plot type:",
                                                             c("-" = "-")
                                                 ),
                                                 radioButtons("findMarkersLabels", label = "Show cluster labels: ",
                                                              choices = list("Yes" = TRUE, 
                                                                             "No" = FALSE)
                                                 ),
                                                 radioButtons("findMarkersOrder", label = "Prioritize expressing cells: ",
                                                              choices = list("Yes" = TRUE, 
                                                                             "No" = FALSE)
                                                 ),
                                                 sliderInput("findMarkersMaxCutoff", "Set max expression value: (quantile)", min = 0, max = 99, value = 99, step = 1),
                                                 sliderInput("findMarkersMinCutoff", "Set minimum expression value: (quantile)", min = 0, max = 99, value = 0, step = 1),
                                                 actionButton(inputId = "findMarkersFPConfirm", label = "Display plot", class="btn-run", icon = icon("check-circle")),
                                                 tags$hr(),
                                                 tags$h3("Add a new signature"),
                                                 textInput(inputId = "findMarkersSignatureName", label = "Gene signature name :", value = "Signature1"),
                                                 textAreaInput(inputId = "findMarkersSignatureMembers", label = "Paste a list of genes", cols = 80, rows = 15, placeholder = "Prg4\nTspan15\nCol22a1\nHtra4"),
                                                 actionButton(inputId = "findMarkersSignatureAdd", label = "Calculate signature score!", class="btn-run", icon = icon("check-circle"))
                                             ),
                                             box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                 div(id="findMarkersFeaturePlot_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "findMarkersFeaturePlot", height = "1300px")
                                                     )
                                                 )
                                             )
                                           )),
                                           tabPanel("Multi-feature vizualization", fluidRow(
                                             box(width=3, status="info", solidHeader=T, title="Options",
                                                 selectizeInput(inputId = 'findMarkersFeaturePair1',
                                                                label = 'Select 1st feature:',
                                                                choices = NULL,
                                                                selected = NULL,
                                                                multiple = FALSE),
                                                 selectizeInput(inputId = 'findMarkersFeaturePair2',
                                                                label = 'Select 2nd Feature:',
                                                                choices = NULL,
                                                                selected = NULL,
                                                                multiple = FALSE),
                                                 sliderInput("findMarkersBlendThreshold", "Select threshold for blending:", min = 0, max = 1, value = 0.5, step = 0.1),
                                                 selectInput("findMarkersFeaturePairReductionType", "Plot type:",
                                                             c("-" = "-")
                                                 ),
                                                 radioButtons("findMarkersFeaturePairLabels", label = "Show cluster labels: ",
                                                              choices = list("Yes" = TRUE, 
                                                                             "No" = FALSE)
                                                 ),
                                                 radioButtons("findMarkersFeaturePairOrder", label = "Prioritize expressing cells: ",
                                                              choices = list("Yes" = TRUE, 
                                                                             "No" = FALSE)
                                                 ),
                                                 sliderInput("findMarkersFeaturePairMaxCutoff", "Set max expression value: (quantile)", min = 0, max = 99, value = 99, step = 1),
                                                 sliderInput("findMarkersFeaturePairMinCutoff", "Set minimum expression value: (quantile)", min = 0, max = 99, value = 0, step = 1),
                                                 actionButton(inputId = "findMarkersFeaturePairConfirm", label = "Display plot", class="btn-run", icon = icon("check-circle"))
                                             ),
                                             box(width=9, status="info", solidHeader=TRUE, title="Plot",
                                                 div(
                                                   column(
                                                     div(id="findMarkersFPfeature1_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId="findMarkersFPfeature1", height = "650px")
                                                         )
                                                     ), width = 6),
                                                   column(
                                                     div(id="findMarkersFPfeature2_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId="findMarkersFPfeature2", height = "650px")
                                                         )
                                                     ), width = 6),
                                                   column(
                                                     div(id="findMarkersFPfeature1_2_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId="findMarkersFPfeature1_2", height = "650px")
                                                         )
                                                     ), width = 6),
                                                   column(
                                                     div(id="findMarkersFPcolorbox_loader",
                                                         shinycssloaders::withSpinner(
                                                           plotlyOutput(outputId="findMarkersFPcolorbox", height = "650px")
                                                         )
                                                     ), width = 6),
                                                 )
                                             )
                                           )
                                           ),
                                           tabPanel("Violin plot", fluidRow(
                                             box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                 selectizeInput(inputId = 'findMarkersGeneSelect2',
                                                                label = 'Type the name of a gene, gene signature or numeric metadata column: 
                                                                      (e.g. "Ccl2" or "Signature1_Ucell", or "nCount_RNA")',
                                                                choices = NULL,
                                                                selected = NULL,
                                                                multiple = FALSE), # allow for multiple inputs
                                                 actionButton(inputId = "findMarkersViolinConfirm", label = "Display plot", class="btn-run", icon = icon("check-circle"))
                                             ),
                                             box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                 div(id="findMarkersViolinPlot_loader",
                                                     shinycssloaders::withSpinner(
                                                       plotlyOutput(outputId = "findMarkersViolinPlot", height = "800px")
                                                     )
                                                 )
                                             )
                                           )
                                           )
                               )
                             )
                    )
        )
)
#tab
    ))
  )
,
server = function(input, output, session){
  
  session$sendCustomMessage("handler_disableTabs", "sidebarMenu") # disable all tab panels (except Data Input) until files are uploaded
  session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "utilities")) #enable GRN tab only for RNA
  options(shiny.maxRequestSize=100*1024^3) 
  hideTab(inputId="grnTabPanel", target="scATAC-seq", session = session)
  
  hideAllLoaders() # helper function (in global.R) that initially hides all loaders
  metaD <- reactiveValues(my_project_name="-", all_lin="0")
  
  
  #Upload tab
     observeEvent(input$upload10xRNAConfirm, {
    if(!is.null(proj_default) | !is.null(seurat_object))
    {
      showModal(modal_confirm)
    }
     else
     {
      session$sendCustomMessage("handler_disableTabs", "sidebarMenu")
      session$sendCustomMessage("handler_disableAllButtons", "upload10xRNAConfirm")
       #tryCatch({
       # showModal(modalDialog(div('Data upload in progress. Please wait...'))) 
         showModal(modalDialog(div('Data upload in progress. Please wait...')))
         metaD$my_project_name <- input$upload10xRNAprojectID
        minimum_cells <<- input$upload10xRNAminCells
        minimum_features <<- input$upload10xRNAminFeatures
        organism <<- input$upload10xRNARadioSpecies

        userId <- session$token
        user_dir <<- paste0("./usr_temp/", userId, metaD$my_project_name, gsub(pattern = "[ ]|[:]", replacement = "_", x = paste0("_", Sys.time())))
        dir.create(user_dir)
        file.copy(from = input$matrix$datapath, to = paste0(user_dir, "/matrix.mtx.gz"), overwrite = TRUE) #,input$matrix$name
        file.copy(from = input$barcodes$datapath, to = paste0(user_dir, "/barcodes.tsv.gz"), overwrite = TRUE) #, input$barcodes$name
        file.copy(from = input$genes$datapath, to = paste0(user_dir, "/features.tsv.gz"), overwrite = TRUE) #, input$genes$name

        seurat_data <- Read10X(user_dir)
        seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))

        init_seurat_object <<- seurat_object

        if(organism == "mouse")
        {
          seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^mt-")
          init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^mt-")
        }
        else
        {
          seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^MT-")
          init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^MT-")
        }

        updateMetadata()
        #updateInputLRclusters1()
        updateSelInpColor()
        updateGeneSearchFP()
        updateQC_choices()
        updateMaxHVGs()
        updatePCs_selection(10)
        updateUtilitiesAssays()
        #disableTabsATAC()
        #update signature text area
        if(organism == "human")
          updateTextAreaInput(session, inputId="findMarkersSignatureMembers", placeholder = "PRG4\nTSPAN15\nCOL22A1\nHTRA4")

        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "2. Droplet processing", "Utilities"))
     #  }, error = function(e) {
     #    print(paste("Error :  ", e))
     #    session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
     #  }, finally = { # with or without error
     #    removeModal()
     #    Sys.sleep(1) # giving some time for rendering for smoother transition
     #    session$sendCustomMessage("handler_enableAllButtons", "upload10xRNAConfirm")
     # })
     }
   })
   
     observeEvent(input$uploadh5Confirm, {
       if(!is.null(proj_default) | !is.null(seurat_object))
       {
         showModal(modal_confirm)
       }
       else
       {
         session$sendCustomMessage("handler_disableTabs", "sidebarMenu") 
         session$sendCustomMessage("handler_disableAllButtons", "uploadh5Confirm") 
         tryCatch({
           showModal(modalDialog(div('Data upload in progress. Please wait...'))) #position:absolute;top:50%;left:50%
           # Create the user directory for the input and output of the analysis
           metaD$my_project_name <- input$uploadh5projectID
           minimum_cells <<- input$uploadh5minCells
           minimum_features <<- input$uploadh5minFeatures
           organism <<- input$uploadh5RadioSpecies
           
           testMatrix <- Read10X_h5(input$uploadh5$datapath)
           seurat_object <<- CreateSeuratObject(counts = testMatrix, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
           init_seurat_object <<- seurat_object
           
           if(organism == "mouse")
           {
             seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^mt-")
             init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^mt-")
           }
           else
           {
             seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^MT-")
             init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^MT-")
           }
           
           updateMetadata()
           #updateInputLRclusters1()
           updateGeneSearchFP()
           updateQC_choices()
           updateMaxHVGs()
           updatePCs_selection(10)
           updateUtilitiesAssays()
     
           #disableTabsATAC()
           #update signature text area
           if(organism == "human")
             updateTextAreaInput(session, inputId="findMarkersSignatureMembers", placeholder = "PRG4\nTSPAN15\nCOL22A1\nHTRA4")
           
           session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "2. Droplet processing","Utilities"))
           updateSelInpColor()
         }, error = function(e) {
           print(paste("Error :  ", e))
           session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
         }, finally = { # with or without error
           removeModal()
           Sys.sleep(1) # giving some time for rendering, for smoother transition
           session$sendCustomMessage("handler_enableAllButtons", "uploadh5Confirm")
         })
       }
     })
   observeEvent(input$uploadSeuratRdsConfirm, {
     if(!is.null(proj_default) | !is.null(seurat_object))
     {
       showModal(modal_confirm)
     }
     else
     {
       session$sendCustomMessage("handler_disableAllButtons", "uploadSeuratRdsConfirm")
     }
     tryCatch({
       showModal(modalDialog(div('Data upload in progress. Please wait...')))
       seurat_object <<- readRDS(input$uploadRdsFile$datapath)
       init_seurat_object <<- seurat_object 
       organism <<- input$uploadRdsRadioSpecies
       print(seurat_object)
       print(init_seurat_object)
       if(is.null(seurat_object@meta.data)) {
         print("Metadata table is missing")
       } else {
         #disableTabsATAC()
         updateMetadata()
         #updateInputLRclusters1()
         updateGeneSearchFP()
         updateQC_choices()
         updateMaxHVGs()
         updatePCs_selection(10)
         updateUtilitiesAssays()
         updateRegressOut()
         #CHECK percent.mt
         if(is.null(seurat_object@meta.data$percent.mt))
         {
           #calculate percent.mt and enable 
           if(organism == "mouse")
           {
             seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^mt-")
             init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^mt-")
           }
           else
           {
             seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^MT-")
             init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^MT-")
           }
         }
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "2. Droplet processing", "Utilities"))
         updateSelInpColor()
         tableMeta <- seurat_object@meta.data
         f <- sapply(tableMeta, is.factor)
         factors <- colnames(tableMeta[, f])
         updateSelectInput(session, "utilitiesActiveClusters", choices = factors)
       }
       
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
     }, finally = { # with or without error
       removeModal()
       Sys.sleep(1) # giving some time for rendering for smoother transition
       session$sendCustomMessage("handler_enableAllButtons", "uploadSeuratRdsConfirm")
     })
   })
   #------------------Doublet detection------------------------------------------
   
   observeEvent(input$doubletsConfirm, {
     
     start.time <- Sys.time()
     print(paste0("!!!!-----Doublets_start-----", start.time))
     session$sendCustomMessage("handler_disableAllButtons", "doubletsConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else
       {
         showModal(modalDialog(div('Doublet detection in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         # if(is.null(input$samplecolumn)){
         #   seurat_object <<- scDblFinder(as.SingleCellExperiment(seurat_object))
         # } else{
         # if(is.null(input$doubletsart)){
         # seurat_object <<- scDblFinder(as.SingleCellExperiment(seurat_object), dbr= input$doubletsPN, nfeatures = 1352, samples=input$samplecolumn)
         # }
         # else{
         #   seurat_object <<- scDblFinder(as.SingleCellExperiment(seurat_object), dbr= input$doubletsPN, nfeatures = input$doubletsart, samples=input$samplecolumn)
         # }
         # }
           if(input$samplecolumn==" "){
         seu<<-scDblFinder(as.SingleCellExperiment(seurat_object))
         } else{
           seu <<- scDblFinder(as.SingleCellExperiment(seurat_object), samples=input$samplecolumn)
         }
         logcounts(seu) <<- assay(seu, "counts")
         seurat_object <<- as.Seurat(seu)
         init_seurat_object <<- seurat_object
         df_total <- table(seurat_object$scDblFinder.class)
         print(head(seurat_object@meta.data))
         updateMetadata()
         updateSelInpColor()
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "3. Quality control", "Utilities"))
         output$doubletsInfo <- renderPrint(
           {
             cat(paste0("Total number of cells: ", nrow(seurat_object@meta.data), "\nTotal number of doublets: ", df_total["doublet"], "\nTotal number of singlets: ", df_total["singlet"]))
           })
       }
     }, error = function(e) 
     {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error during doublet detection.")
     }, finally = 
       {
         end.time <- Sys.time()
         time.taken <- as.numeric (end.time - start.time, units = "mins")
         print(paste0("!!!!-----Doublets_end_Mins_taken-----", time.taken))
         
         removeModal()
         Sys.sleep(1)
         session$sendCustomMessage("handler_enableAllButtons", "doubletsConfirm")
       })
     
     
   })
   
   observeEvent(input$doubletsRemove, { 
     session$sendCustomMessage("handler_disableAllButtons", "doubletsRemove")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else
       {
         seurat_object <<- subset(seurat_object, subset = scDblFinder.class=="singlet")
         init_seurat_object<<- seurat_object
         updateMetadata()
         shinyalert("Doublets were removed successfully.", type = "success")
       }
     }, error = function(e) 
     {
       print(paste("Error :  ", e))
       shinyalert("There was an error during doublets removal. Maybe doublet detection is not completed?", type = "error")
     }, finally = 
       {
         session$sendCustomMessage("handler_enableAllButtons", "doubletsRemove")
       }
     )
   })
   ##QC tab
   
   observeEvent(input$qcDisplay, {
     start.time <- Sys.time()
     print(paste0("!!!!-----DisplayQC_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "qcDisplay")
     tryCatch({
       if (identical(init_seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         shinyjs::show("nFeatureViolin_loader")
         shinyjs::show("totalCountsViolin_loader")
         shinyjs::show("mitoViolin_loader")
         shinyjs::show("genesCounts_loader")
         shinyjs::show("mtCounts_loader")
         showModal(modalDialog(div('Quality control in progress. This operation may take a few minutes for large datasets. Please wait...')
                               )) #position:absolute;top:50%;left:50%

         output$nFeatureViolin <- renderPlotly(
           {
             p <- VlnPlot(init_seurat_object, features = c("nFeature_RNA"), pt.size = 0, group.by = "orig.ident",
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) +
               theme_bw() +
               geom_hline(yintercept=c(as.numeric(input$minUniqueGenes), as.numeric(input$maxUniqueGenes)), linetype="dashed", color = "red", linewidth=1) +
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none") +
               labs(title = "", y="Genes detected/cell")

             p_temp <- VlnPlot(init_seurat_object, features = c("nFeature_RNA"), pt.size = 0, group.by = "orig.ident",
                               cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) +
               theme_bw() +
               geom_hline(yintercept=c(as.numeric(input$minUniqueGenes), as.numeric(input$maxUniqueGenes)), linetype="dashed", color = "red", linewidth=1) +
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none") +
               labs(title = "", y="Detected genes/cell")

             plotly::ggplotly(p + ylim(min(init_seurat_object$nFeature_RNA)-1, max(init_seurat_object$nFeature_RNA)+2), tooltip = c("x", "y"))
           }
         )

         output$totalCountsViolin <- renderPlotly(
           {
             p <- VlnPlot(init_seurat_object, features = c("nCount_RNA"), pt.size = 0, group.by = "orig.ident",
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) +
               theme_bw() +
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none")+
               labs(title = "", y="Total counts/cell")
             plotly::ggplotly(p, tooltip = c("x", "y"))
           }
         )

         output$mitoViolin <- renderPlotly(
           {
             p <- VlnPlot(init_seurat_object, features = c("percent.mt"), pt.size = 0, group.by = "orig.ident",
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) +
               theme_bw() +
               #scale_y_continuous(labels = number_format(accuracy = 0.001)) +
               geom_hline(yintercept= as.numeric(input$maxMtReads), linetype="dashed", color = "red", linewidth=1) +
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none")+
               labs(title = "", y="% of reads mapped to mitochondrial genome/cell")
             plotly::ggplotly(p, tooltip = c("x", "y")) #%>%
             #layout(yaxis = list(hoverformat = ".2f"))
           }
         )

         output$mtCounts <- renderPlotly(
           {
             p <- FeatureScatter(init_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", raster = F,
                                 cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident']))))
             gp <- plotly::ggplotly(p)
             print(gp)
           }
         )

         output$genesCounts <- renderPlotly(
           {
             p <- FeatureScatter(init_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", raster = F,
                                 cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident']))))
             gp <- plotly::ggplotly(p)
             print(gp)
           }
         )

         output$cellStats <- renderPrint(
           {
             cat(paste0("Total number of cells: ", nrow(init_seurat_object@meta.data)))
           }
         )
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "The selected Quality Control arguments cannot produce meaningful visualizations.")
     }, finally = {
       removeModal()

       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----DisplayQC_end_Mins_taken-----", time.taken))

       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "qcDisplay")
     })
   })
   
   observeEvent(input$qcConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----QCFiltering_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "qcConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         updateTabsetPanel(session, inputId = "qc_tabs_rna", selected = "Post-filtering plots")
         showModal(modalDialog(div('Analysis in Progress. This operation may take several minutes, please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("filteredNFeatureViolin_loader")
         shinyjs::show("filteredTotalCountsViolin_loader")
         shinyjs::show("filteredMitoViolin_loader")
         shinyjs::show("filteredGenesCounts_loader")
         shinyjs::show("filteredMtCounts_loader")
         
         qc_minFeatures <<- input$minUniqueGenes
         qc_maxFeatures <<- input$maxUniqueGenes
         qc_maxMtPercent <<- input$maxMtReads
         
         seurat_object <<- subset(init_seurat_object, subset = nFeature_RNA > as.numeric(qc_minFeatures) & nFeature_RNA < as.numeric(qc_maxFeatures) & percent.mt < as.double(qc_maxMtPercent)) #filter object
         
         output$filteredNFeatureViolin <- renderPlotly(
           {
             p <- VlnPlot(seurat_object, features = c("nFeature_RNA"), pt.size = 0, group.by = input$qcColorBy,
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
               theme_bw() + 
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none") +
               labs(title = "", y="Genes detected/cell")
             plotly::ggplotly(p, tooltip = c("x", "y")) 
           }
         )
         output$filteredTotalCountsViolin <- renderPlotly(
           {
             p <- VlnPlot(seurat_object, features = c("nCount_RNA"), pt.size = 0, group.by = input$qcColorBy,
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
               theme_bw() + 
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none")+
               labs(title = "", y="Total counts/cell")
             plotly::ggplotly(p, tooltip = c("x", "y")) 
           }
         )
         output$filteredMitoViolin <- renderPlotly(
           {
             p <- VlnPlot(seurat_object, features = c("percent.mt"), pt.size = 0, group.by = input$qcColorBy,
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
               theme_bw() + 
               theme(
                 plot.title = element_blank(),
                 axis.title.x = element_blank(),
                 legend.position = "none")+
               labs(title = "", y="% of reads mapped to mitochondrial genome/cell")
             plotly::ggplotly(p, tooltip = c("x", "y"))
           }
         )
         output$filteredMtCounts <- renderPlotly(
           {
             p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = input$qcColorBy, raster = F, 
                                 cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
             gp <- plotly::ggplotly(p)
             print(gp)
           }
         )
         output$filteredGenesCounts <- renderPlotly(
           {
             p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = input$qcColorBy, raster = F,
                                 cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
             gp <- plotly::ggplotly(p)
             print(gp)
           }
         )
         output$filteredCellStats <- renderPrint(
           {
             cat(paste0("\nTotal number of cells after filtering: ", nrow(seurat_object@meta.data)))
           }
         )
         updateSelInpColor()
         session$sendCustomMessage("handler_disableTabs", "sidebarMenu")
         updateMetadata()
         cleanAllPlots(F) # fromDataInput -> FALSE
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "2. Droplet processing","3. Quality control", "4. Data normalization", "Utilities"))
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "The selected Quality Control arguments cannot produce meaningful visualizations.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----QCfilter_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "qcConfirm")
     })
   })
  
   
   #------------------Normalization tab--------------------------------
   observeEvent(input$normalizeConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----Scaling_start-----", start.time))
     session$sendCustomMessage("handler_log", "Starting normalization procedure")
     session$sendCustomMessage("handler_disableAllButtons", "normalizeConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         shinyjs::show("hvgScatter_loader")
         showModal(modalDialog(div('sctransform normalization in progress. This operation can take a while.
                                  Please wait...')))
         normalize_scaleRegressOut <- input$normalizeRegressColumns
         if(is.null(normalize_scaleRegressOut)) seurat_object <<- SCTransform(seurat_object, vst.flavor = "v2")
         else seurat_object <<- SCTransform(seurat_object, vst.flavor = "v2", vars.to.regress= normalize_scaleRegressOut)

         # Rendering
         output$hvgScatter <- renderPlotly({ # tooltip
           if(length(VariableFeatures(seurat_object)) != 0){
             plot1 <- VariableFeaturePlot(seurat_object)
             varplot <- plot1$data
             varplot$gene <- rownames(varplot)
             varplot$colors[varplot$colors == "yes"] <- paste0("Highly Variable genes(", length(VariableFeatures(seurat_object)), ")")
             varplot$colors[varplot$colors == "no"] <- paste0("Not Variable genes(", length(rownames(seurat_object)) - length(VariableFeatures(seurat_object)), ")")
             p <- ggplot(varplot, aes(x=log10(gmean), y=residual_variance, color=colors, text = paste0("log10(Geometric Mean of Expression): ", log10(gmean),
                                                                                                       "\nResidual Variance: ", residual_variance,
                                                                                                       "\ngene: ", gene)
             )) +
               geom_point()+
               theme_bw() +
               scale_color_manual(
                 values = c("red", "black")
               )+
               labs(x="log10(Geometric Mean of Expression)", y="Residual Variance", color="")
             gp <- plotly::ggplotly(p, tooltip = "text")
             return(gp)
           } else session$sendCustomMessage("handler_alert", "There are no variable features in the Seurat object.")
         }) # End Rendering
         output$hvgTop10Stats <- renderPrint(
           {
             top10hvgs <- head(VariableFeatures(seurat_object), 10)
             cat("\nTop10 variable features:", c(top10hvgs[1:10]))
           })
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "5. PCA"))
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "The selected Normalization arguments cannot produce meaningful visualizations.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----Scaling_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       print("testing")
       session$sendCustomMessage("handler_enableAllButtons", "normalizeConfirm")
       
       })
   })
   
   #------------------PCA tab------------------------------------------
   observeEvent(input$PCrunPCA, {
     start.time <- Sys.time()
     print(paste0("!!!!-----PCA_start-----", start.time))
     session$sendCustomMessage("handler_disableAllButtons", "PCrunPCA")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!("ScaleData.RNA" %in% names(seurat_object@commands)|"SCTransform.RNA"%in% names(seurat_object@commands))) session$sendCustomMessage("handler_alert", "Data need to be normalized and scaled first. Please, execute the normalization step.")
       else {
         showModal(modalDialog(div('PCA analysis in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("elbowPlotPCA_loader")
         shinyjs::show("PCAscatter_loader")
         seurat_object <<- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 50)
          updateUmapTypeChoices("pca")
          
          output$elbowPlotPCA <- renderPlotly(
            {
              plot1 <- ElbowPlot(seurat_object, ndims = 50)
              plot1_data <- plot1$data
              colnames(plot1_data)[1] <- "PC"
              colnames(plot1_data)[2] <- "SD"
              p <- ggplot(plot1_data, aes(x=PC, y=SD))+
                geom_point() +
                theme_bw() +
                labs(x="PC", y="Standard Deviation")
              gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
              print(gp)
            })
          
          output$PCAscatter <- renderPlotly(
            {
              #prepare metadata
              meta <- seurat_object@meta.data
              meta$Cell_id <- rownames(meta)
              meta <- meta[, ]#meta[, c('Cell_id', 'seurat_clusters', 'orig.ident')]
              reduc_data <- data.frame()
              #prepare colors
              cols = colorRampPalette(brewer.pal(12, "Paired"))(1)
              #umap data frame
              seurat_object_reduc <- as.data.frame(seurat_object@reductions$pca@cell.embeddings)
              seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
              seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
              reduc_data <- left_join(seurat_object_reduc, meta)
              
              p <- ggplot(data=reduc_data, aes_string(x="PC_1", y="PC_2", fill="orig.ident")) +
                geom_point(size=2, shape=19, stroke=0)+
                scale_fill_manual(values = cols)+
                scale_size()+
                theme_bw() +
                theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                      axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                      axis.title.y = element_text(face = "bold", color = "black", size = 25),
                      axis.title.x = element_text(face = "bold", color = "black", size = 25),
                      legend.text = element_text(face = "bold", color = "black", size = 9),
                      legend.title = element_text(face = "bold", color = "black", size = 9),
                      legend.position="right",
                      title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
                labs(x="PC 1", y="PC 2", title = "", fill="Color")
              print(plotly::ggplotly(p))
            })
          session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "6. Clustering"))
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the PCA analysis.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----PCA_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "PCrunPCA")
     }
     )
   })
   observeEvent(input$PCconfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "PCconfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first Run PCA above.")
       else {
         shinyjs::show("PCAloadings_loader")
         shinyjs::show("PCAheatmap_loader")
         showModal(modalDialog(div('Rendering in progress, please wait...')) )

         activePC <- as.numeric(input$PCin)
       #   
         output$PCAloadings <- renderPlotly(
           {
             plot1 <- VizDimLoadings(seurat_object, dims = activePC, reduction = "pca", balanced = TRUE)
             plot1_data <- plot1$data
             plot1_data$orig.ident <- unique(seurat_object@meta.data$orig.ident)
             colnames(plot1_data)[1] <- "LoadingScore"
             activePC <- paste0("PC", "_", activePC)

             export_loadingScoresTable_RNA <<- plot1$data

             p <- ggplot(plot1_data, aes(x=LoadingScore, y=feature, color=orig.ident))+
               geom_point() +
               theme_bw() +
               scale_color_manual(values="blue")+
               theme(legend.position = "none")+
               labs(x=activePC, y="")

             gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
             print(gp)
           }
         )

         output$PCAheatmap <- renderPlotly(
           {
             p <- DimHeatmap(seurat_object, dims = activePC, cells = 500, balanced = TRUE, fast = F)
             p
             plotly::ggplotly(p)
           }
         )
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the visualization of PCA analysis.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "PCconfirm")
     })
   })
   
   output$pcaRNAExport <- downloadHandler(
     filename = function() { 
       paste("loadingScoresTableRNA-", Sys.Date(), ".txt", sep="")
     },
     content = function(file) {
       write.table(export_loadingScoresTable_RNA, file, sep = "\t", quote = F, row.names = F)
     })
   
   #------------------Clustering tab------------------------------------------
   observeEvent(input$snnConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----Clustering_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "snnConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Perform PCA first, please.")
       else {
         showModal(modalDialog(div('Clustering cells in progress. Please wait...')))
         snn_dims <<- input$snnPCs
         snn_k <<- input$snnK
         cluster_res <<- input$clusterRes
         cluster_dims <<- input$clusterPCs
         
         seurat_object <<- FindNeighbors(seurat_object, k.param = as.numeric(snn_k), dims = 1:as.numeric(snn_dims), reduction = "pca", nn.method="annoy")
         seurat_object <<- FindClusters(seurat_object, resolution = as.numeric(cluster_res))
         
         updateClusterTab()
         updateSelInpColor()
         updateInputLRclusters()
         updateMetadata()
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "7. UMAP"))
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the clustering procedure.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----Clustering_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "snnConfirm")
     })
   })
   
   observeEvent(input$clusterBarplotConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "clusterBarplotConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute clustering first.")
       else {
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("clusterBarplot_loader")
         
         if(input$clusterGroupBy == "seurat_clusters")
         {
           #barplot for cell distribution per cluster
           clusterTable <- as.data.frame(table(seurat_object$seurat_clusters)) #as.data.frame(table(Idents(seurat_object)))
           totalCells <- sum(clusterTable$Freq)
           
           clusterTable$Perc <- (clusterTable$Freq)/totalCells
           colnames(clusterTable)[1] <- "Cluster"
           colnames(clusterTable)[2] <-  "Cells"
           colnames(clusterTable)[3] <- "Percentage"
           cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(clusterTable$Cluster)))
           
           p <- ggplot(clusterTable) + theme_bw() +
             geom_bar( mapping = aes(x = Cluster, y = Percentage, fill=Cluster), stat = "identity" ) +
             scale_fill_manual(values = cols ) +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 12),
                   axis.title.x = element_text(face = "bold", color = "black", size = 12),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black")) +
             labs(x="", y="Percent of cells", fill="Clusters")
           gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
           
           output$clusterBarplot <- renderPlotly({print(gp)})
         }
         else
         {
           cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object$seurat_clusters)))
           
           p <- dittoBarPlot(seurat_object, var = "seurat_clusters",group.by = input$clusterGroupBy, scale = "percent") + theme_bw() +
             scale_fill_manual(values = cols ) +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 12),
                   axis.title.x = element_text(face = "bold", color = "black", size = 12),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black")) 
           
           gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
           
           output$clusterBarplot <- renderPlotly({print(gp)})
         }
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error during plotting the barplot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "clusterBarplotConfirm")
     })
   })
   
   observeEvent(input$snnDisplayConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----SNN_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "snnDisplayConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('SNN graph visualization in progress. This a slow operation for large datasets. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("snnSNN_loader")
         
         output$snnSNN <- renderVisNetwork(
           {
             if(!is.null(seurat_object@graphs$SCT_snn))
             {
               #set.seed(9)
               
               mygraph <- as.matrix(seurat_object@graphs$SCT_snn)
               graphOut <- graph_from_adjacency_matrix(mygraph, mode = "undirected", weighted = T, diag = F)
               graphSimple <- simplify(graphOut) #, remove.loops=T)
               weights <- E(graphSimple)$weight
               sub_nodes <- V(graphSimple)$name
               
               tableCl <- seurat_object@meta.data[, ]
               tableCl$Cell_id <- rownames(tableCl)
               tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
               tableCl <- tableCl[order(tableCl$seurat_clusters), ]
               tableCl <- tableCl[sub_nodes, ]
               colors_cl <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
               tableCl$color <- colors_cl[as.numeric(tableCl$seurat_clusters)]
               
               V(graphSimple)$color <- tableCl$color
               visIgraph(graphSimple, layout = "layout_with_lgl", randomSeed = 9) %>%
                 visInteraction(navigationButtons = TRUE, hover = TRUE)
             }
           }
         )
         
         
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was a problem in drawing the SNN graph.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----SNN_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "snnDisplayConfirm")
     })
   }) 
   
   output$clusterTableRNAExport <- downloadHandler(
     filename = function() { 
       paste("clusterTableRNA-", Sys.Date(), ".txt", sep="")
     },
     content = function(file) {
       write.table(export_clustertable_RNA, file, sep = "\t", quote = F, row.names = F)
     })
   #find markers -----------
   
   observeEvent(input$findMarkersConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----DEA_start-----", start.time))
     # 
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute clustering first .")
       else {
         showModal(modalDialog(div('Finding markers. This operation may take several minutes, please wait...')))
         markers_logFCBase <<- input$findMarkersLogBase

         if(markers_logFCBase == "avg_logFC")
         {
           seurat_object@misc$markers <<- FindAllMarkers(seurat_object, test.use = input$findMarkersTest, min.pct = as.numeric(input$findMarkersMinPct), logfc.threshold = as.numeric(input$findMarkersLogFC),
                                                         return.thresh = as.numeric(input$findMarkersPval), base = exp(1))
         }
         else
         {
           seurat_object@misc$markers <<- FindAllMarkers(seurat_object, test.use = input$findMarkersTest, min.pct = as.numeric(input$findMarkersMinPct), logfc.threshold = as.numeric(input$findMarkersLogFC),
                                                         return.thresh = as.numeric(input$findMarkersPval), base = 2)
         }
         updateInputLRclusters()

         output$findMarkersTable <- renderDataTable(
           {
             if (!is.null(seurat_object@misc$markers))
             {
               seurat_object@misc$markers <<- seurat_object@misc$markers[, c(7, 6, 1:5)]
               output$findMarkersTable <- renderDataTable(seurat_object@misc$markers, options = list(pageLength = 20))
               export_markerGenes_RNA <<- seurat_object@misc$markers
             }
           }
         )

         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "11. Cell-type annotation"))
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the DE Analysis")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("DEA_end_Mins_taken", time.taken))

       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersConfirm")
     })
   })
   
   observeEvent(input$findMarkersTop10HeatmapConfirm, {
     #DEA output rendering
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersTop10HeatmapConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersHeatmap_loader")
         
         top10 <- seurat_object@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = eval(parse(text=markers_logFCBase))) 
         
         downsampled <- seurat_object
         if(ncol(seurat_object) > 1500)
         {
           set.seed(9)
           downsampled <- subset(seurat_object, cells = sample(Cells(seurat_object), 1500))  
         }
         
         scaled_tabe <- as.data.frame(downsampled@assays$SCT@scale.data)
         scaled_tabe$gene <- rownames(scaled_tabe)
         scaled_tabe_order <- as.data.frame(top10$gene)
         colnames(scaled_tabe_order)[1] <- "gene"
         scaled_tabe_final <- left_join(scaled_tabe_order, scaled_tabe)
         scaled_tabe_final <- na.omit(scaled_tabe_final)
         tableCl <- downsampled@meta.data[, ]
         tableCl$Cell_id <- rownames(tableCl)
         tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
         tableCl <- tableCl[order(tableCl$seurat_clusters), ]
         
         clip<-function(x, min=-2, max=2) {
           x[x<min]<-min; 
           x[x>max]<-max; 
           x
         }
         
         final_mat <- scaled_tabe_final[, -1]
         final_mat <- final_mat[, tableCl$Cell_id]
         
         cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
         names(cols) <- unique(tableCl$seurat_clusters)
         
         output$findMarkersHeatmap <- renderPlotly({
           heatmaply(clip(final_mat), Rowv = F, Colv = F, colors = rev(RdBu(256)), showticklabels = c(F, T), labRow  = scaled_tabe_final$gene, 
                     col_side_colors = tableCl$seurat_clusters, col_side_palette =  cols)
         })
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the visualization of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersTop10HeatmapConfirm")
     })  
     
   })
   
   observeEvent(input$findMarkersTop10DotplotConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersTop10DotplotConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data  first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersDotplot_loader")
         output$findMarkersDotplot <- renderPlotly(
           {
             top10 <- seurat_object@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = eval(parse(text=markers_logFCBase))) 
             p <- DotPlot(seurat_object, features = rev(unique(top10$gene)), dot.scale = 6, cols = c("grey", "red")) + RotatedAxis() + labs(fill="Average\nexpression")
             plotly::ggplotly(p)
           }
         )
         
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the visualization of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersTop10DotplotConfirm")
     })
   })
   
   observeEvent(input$findMarkersSignatureAdd, {
     start.time <- Sys.time()
     print(paste0("!!!!-----Signature_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersSignatureAdd")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Gene signature addition in progress. Please wait...')))
         markers <- list()
         varTextarea <- input$findMarkersSignatureMembers
         markers2 <- as.vector(unlist(strsplit(varTextarea, "\\n"))) 
         print(markers2)
         print(class(markers2))
         print(length(markers2))
         markers[[1]] <- unlist(strsplit(varTextarea, "\\n"))
         sig_name <- input$findMarkersSignatureName
         print(sig_name)
         names(markers)[1] <- sig_name
         
         genesFound <- markers2[(which(markers2 %in% rownames(seurat_object)))]
         print(genesFound)
         
         if(length(genesFound) == 0)
         {
           print("inside all zero")
           session$sendCustomMessage("handler_alert", "None of the selected genes were found in the dataset.The signature scores could not be calculated.")
         }
         else
         {
           seurat_object <<- AddModuleScore_UCell(seurat_object, features = markers)
           updateGeneSearchFP()
           updateMetadata()
           if(length(genesFound) != length(markers2))
           {
             print("inside some not found")
             session$sendCustomMessage("handler_alert", "Some of the selected genes were not found in the dataset and are excluded from the signature's calculation.")
           }
         }
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "The signature could not be added.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----Signature_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersSignatureAdd")
     })
   })
   observeEvent(input$findMarkersFPConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersFPConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersFeaturePlot_loader")
         
         if(input$findMarkersReductionType != "-")
         {
           geneS <- input$findMarkersGeneSelect
           
           label_x <- ""
           label_y <- ""
           show_label <- as.logical(input$findMarkersLabels)
           order_exp <- as.logical(input$findMarkersOrder)
           minq <- paste0("q", input$findMarkersMinCutoff)
           maxq <- paste0("q", input$findMarkersMaxCutoff)
           
           if(input$findMarkersReductionType == "umap")
           {
             label_x <- "umap_1"
             label_y <- "umap_2"
           }
           else if(input$findMarkersReductionType == "tsne")
           {
             label_x <- "tSNE_1"
             label_y <- "tSNE_2"
           }
           else if(input$findMarkersReductionType == "dfm")
           {
             label_x <- "DC_1"
             label_y <- "DC_2"
           }
           else if(input$findMarkersReductionType == "pca")
           {
             label_x <- "PC_1"
             label_y <- "PC_2"
           }
           else if(input$findMarkersReductionType == "phate")
           {
             label_x <- "PHATE_1"
             label_y <- "PHATE_2"
           }
           
           plot_temp <- FeaturePlot(seurat_object, features = geneS, pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red"), 
                                    order = order_exp, reduction = input$findMarkersReductionType, max.cutoff = maxq, min.cutoff = minq)
           
           plot <- FeaturePlot(seurat_object, features = geneS, pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red"), 
                               order = order_exp, reduction = input$findMarkersReductionType, max.cutoff = maxq, min.cutoff = minq, raster = F) + 
             xlim(min(plot_temp$data[label_x]), max(plot_temp$data[label_x])) + ylim(min(plot_temp$data[label_y]), max(plot_temp$data[label_y])) +
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x=label_x, y=label_y, title = geneS, color="")
           gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
           output$findMarkersFeaturePlot <- renderPlotly({
             gp
           })
         }
         
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was a problem with the generation of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersFPConfirm")
     })
   })
   
   observeEvent(input$findMarkersFeaturePairConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersFeaturePairConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersFPfeature1_loader")
         shinyjs::show("findMarkersFPfeature2_loader")
         shinyjs::show("findMarkersFPfeature1_2_loader")
         shinyjs::show("findMarkersFPcolorbox_loader")
         
         updateFeaturePair()
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was a problem with the generation of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersFeaturePairConfirm")
     })
   })
   
   observeEvent(input$findMarkersViolinConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersViolinConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersViolinPlot_loader")
         if (!identical(input$findMarkersGeneSelect2, NULL)){
           output$findMarkersViolinPlot <- renderPlotly(
             {
               
               geneS <- input$findMarkersGeneSelect2
               print(geneS)
               plot <- VlnPlot(seurat_object, features = geneS, pt.size = 0, raster = F,
                               cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'seurat_clusters'])))) +
                 theme_bw() +
                 theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                       axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                       axis.title.y = element_text(face = "bold", color = "black", size = 25),
                       axis.title.x = element_text(face = "bold", color = "black", size = 25),
                       legend.text = element_text(face = "bold", color = "black", size = 9),
                       legend.title = element_text(face = "bold", color = "black", size = 9),
                       title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
                 labs(x="Cluster", y="", title = geneS, fill="Cluster")
               gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
               gp
             }
           )
         }
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was a problem with the generation of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersViolinConfirm")
     })
   })
   
   observeEvent(input$findMarkersVolcanoConfirm, {
     session$sendCustomMessage("handler_disableAllButtons", "findMarkersVolcanoConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else{
         showModal(modalDialog(div('Rendering in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("findMarkersVolcanoPlot_loader")
         diff_exp_genes <- seurat_object@misc$markers
         cluster_degs <- diff_exp_genes[which(diff_exp_genes$cluster == input$findMarkersClusterSelect), ]
         cluster_degs$status <- "Down regulated"
         for(i in 1:length(cluster_degs$gene))
         {
           if(cluster_degs[i, markers_logFCBase] > 0) 
           {
             cluster_degs$status[i] <- "Up regulated"
           }
         }
         cluster_degs$log10Pval <- -log10(cluster_degs$p_val)
         
         p <- ggplot(data=cluster_degs, aes_string(x=markers_logFCBase, y="log10Pval", fill="status", label="gene", color="status")) + 
           geom_point(size=1, shape=16)+
           scale_fill_manual(values = c("cyan3", "orange"))+
           scale_color_manual(values = c("cyan3", "orange"))+
           scale_size()+
           theme_bw() +
           theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                 axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                 axis.title.y = element_text(face = "bold", color = "black", size = 25),
                 axis.title.x = element_text(face = "bold", color = "black", size = 25),
                 legend.text = element_text(face = "bold", color = "black", size = 9),
                 legend.title = element_text(face = "bold", color = "black", size = 9),
                 legend.position="right",
                 title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
           labs(x=paste0("", markers_logFCBase), y="-log10(Pvalue)", fill="Color", color="")
         
         output$findMarkersVolcanoPlot <- renderPlotly(
           {
             plotly::ggplotly(p, tooltip = c("x", "y", "label"))
           })
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was a problem with the generation of the plot.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "findMarkersVolcanoConfirm")
     })
   })
   
   #umap
   observeEvent(input$umapRunUmap, {
     start.time <- Sys.time()
     print(paste0("!!!!-----UMAP_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "umapRunUmap")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunTsne")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunDFM")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunPhate")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, run PCA first.")
       else {
         showModal(modalDialog(div('Calculating UMAP. Please wait...'))) #position:absolute;top:50%;left:50%
         seurat_object <<- RunUMAP(seurat_object, dims = 1:as.numeric(input$umapPCs), seed.use = as.numeric(input$umapSeed), n.components = 2, reduction = "pca")
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "8. Cell cycle analysis",  "9. Dataset inspection", "10. Marker gene detection"))
         updateUmapTypeChoices("umap")
         updateReduction2("umap")
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the UMAP calculation.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----UMAP_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapRunUmap")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunTsne")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunDFM")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunPhate")
     })
   })
   
   observeEvent(input$umapRunTsne, {
     session$sendCustomMessage("handler_disableAllButtons", "umapRunUmap")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunTsne")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunDFM")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunPhate")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, run PCA first.")
       else {
         showModal(modalDialog(div('tSNE reduction is calculated! Please wait...'))) #position:absolute;top:50%;left:50%
       
           seurat_object <<- RunTSNE(seurat_object, dims = 1:as.numeric(input$umapPCs), seed.use = as.numeric(input$umapSeed), dim.embed = 2, reduction = "pca", verbose = T)
           
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "9. Dataset inspection","10. Marker gene detection"))
         updateUmapTypeChoices("tsne")
         updateReduction2("tsne")
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the tSNE calculation.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapRunUmap")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunTsne")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunDFM")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunPhate")
     })
   })
   
   observeEvent(input$umapRunDFM, {
     session$sendCustomMessage("handler_disableAllButtons", "umapRunUmap")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunTsne")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunDFM")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunPhate")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, run PCA first.")
       else {
         showModal(modalDialog(div('Diffusion map reduction is calculated! Please wait...'))) #position:absolute;top:50%;left:50%
         #prepare input
         dfm_in <- as.data.frame(seurat_object@reductions[['pca']]@cell.embeddings)
         dfm_in$Cell_id <- rownames(dfm_in)
         dfm_in <- dfm_in[, c(51, 1:as.numeric(input$umapPCs))]
         
         #run DFM default
         dfm <- DiffusionMap(dfm_in, n_eigs = 2)
         dfm_out <- dfm@eigenvectors
         colnames(dfm_out) <- gsub("DC", "DC_", colnames(dfm_out))
         
         #add new reduction in seurat_object
         seurat_object[["dfm"]] <<- CreateDimReducObject(embeddings = dfm_out, key = "DC_", assay = DefaultAssay(seurat_object), global = T)
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "9. Dataset inspection","10. Marker gene detection"))
         updateUmapTypeChoices("dfm")
         updateReduction2("dfm")
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the Diffusion Map calculation.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapRunUmap")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunTsne")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunDFM")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunPhate")
     })
   })
   
   observeEvent(input$umapRunPhate, {
     start.time <- Sys.time()
     print(paste0("!!!!-----Phate_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "umapRunUmap")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunTsne")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunDFM")
     session$sendCustomMessage("handler_disableAllButtons", "umapRunPhate")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, run PCA first.")
       else {
         showModal(modalDialog(div('PHATE reduction is calculated! This is slow procedure, please wait...'))) #position:absolute;top:50%;left:50%
         seurat_object <<- RunPHATE(seurat_object, dims = 1:as.numeric(input$umapPCs), n.components = 2, reduction = "pca", seed.use = as.numeric(input$umapSeed))
         session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "9. Dataset inspection","10. Marker gene detection"))
         updateUmapTypeChoices("phate")
         updateReduction2("phate")
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with the Phate calculation.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----Phate_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapRunUmap")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunTsne")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunDFM")
       session$sendCustomMessage("handler_enableAllButtons", "umapRunPhate")
     })
   })
   
   observeEvent(input$umapConfirm, { 
     if(input$umapType != "-")
       updateReduction()
   })
   
   ##Dataset inspection
   
   #------------------Cell cycle tab---------------------------------------------
   observeEvent(input$cellCycleRun, {
     start.time <- Sys.time()
     print(paste0("!!!!-----Cell_cycle_start-----", start.time))

  session$sendCustomMessage("handler_disableAllButtons", "cellCycleRun")
  tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data  first.")
       else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, perform PCA first.")
       else {
         showModal(modalDialog(div('Cell cycle phase analysis in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("cellCyclePCA_loader")
         shinyjs::show("cellCycleBarplot_loader")

         if(organism == "mouse")
         {

           s.genes <- mmus_s
           g2m.genes <- mmus_g2m
         }
         else
         {
           s.genes <- cc.genes.updated.2019$s.genes
           g2m.genes <- cc.genes.updated.2019$g2m.genes
         }

         seurat_object <<- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
         seurat_object@meta.data$CC.Difference <<- seurat_object$S.Score - seurat_object$G2M.Score
         seurat_object@meta.data$Phase <<- factor(seurat_object@meta.data$Phase, levels = c("G1", "S", "G2M"))
         updateMetadata()

         output$cellCyclePCA <- renderPlotly(
           {
             if(!is.null(seurat_object) & input$cellCycleReduction != "-")
             {
               meta <- seurat_object@meta.data
               meta$Cell_id <- rownames(meta)
               label_x <- "PC_1"
               label_y <- "PC_2"
               selected_Reduction <- input$cellCycleReduction

               if(selected_Reduction == "umap")
               {
                 label_x <- "umap_1"
                 label_y <- "umap_2"
               }
               else if(selected_Reduction == "tsne")
               {
                 label_x <- "tSNE_1"
                 label_y <- "tSNE_2"
               }
               else if(selected_Reduction == "dfm")
               {
                 label_x <- "DC_1"
                 label_y <- "DC_2"
               }
               else if(selected_Reduction == "phate")
               {
                 label_x <- "PHATE_1"
                 label_y <- "PHATE_2"
               }

               plot1 <- DimPlot(seurat_object, reduction = selected_Reduction)
               plot1_data <- plot1$data
               plot1_data$Cell_id <- rownames(plot1_data)
               plot1_data <- left_join(plot1_data, meta)

               p <- ggplot(plot1_data, aes_string(x=label_x, y=label_y, color="Phase"))+
                 geom_point() +
                 theme_bw() +
                 labs(x=label_x, y=label_y)+
                 theme(legend.position = "right")

               gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
               print(gp)
             }
           }
         )

         output$cellCycleBarplot <- renderPlotly(
           {
             #barplot
             fsize <- 24

             obj_meta <- seurat_object@meta.data
             obj_meta$cell_id <- row.names(obj_meta)
             obj_meta$Cluster_name <- Idents(seurat_object)
             obj_meta <- obj_meta[, c('cell_id', 'Cluster_name', 'Phase')]
             final_df_seurat <- obj_meta %>% dplyr::group_by(Cluster_name) %>% dplyr::count(Phase)
             colnames(final_df_seurat)[3] <- "Cells"

             total_cells_per_cluster <- as.data.frame(table(seurat_object$seurat_clusters)) #as.data.frame(table(Idents(seurat_object)))
             colnames(total_cells_per_cluster)[1] <- "Cluster_name"
             colnames(total_cells_per_cluster)[2] <- "Total"

             final_df_seurat <- left_join(final_df_seurat, total_cells_per_cluster)
             final_df_seurat$Percentage <- (final_df_seurat$Cells/final_df_seurat$Total)*100

             p <- ggplot(final_df_seurat, aes(x=Phase, y=Percentage, fill=Phase))+
               geom_bar(stat='identity')+
               facet_wrap(~Cluster_name) +
               theme_bw() +
                theme(axis.text.x = element_text(face = "bold", color = "black", size = fsize, angle = 90, vjust = 0.5),
                     axis.text.y = element_text(face = "bold", color = "black", size = 18, angle = 0),
                     axis.title.y = element_text(face = "bold", color = "black", size = fsize),
                     axis.title.x = element_text(face = "bold", color = "black", size = fsize),
                     legend.text = element_text(size = 18, color = "black", face = "bold.italic"),
                     legend.title = element_text(size = 18, color = "black", face = "bold.italic"),
                     strip.text.x = element_text(size = 18, color = "black", face = "bold.italic"),
                       strip.background = element_rect(color="black", size=1.5, linetype="solid")) +
               labs(y="% of cells", x="Phase", color="Phase")
             p
             plotly::ggplotly(p)
           }
         )
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with cell cycle phase analysis.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----Cell_cycle_end_Mins_taken-----", time.taken))

       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "cellCycleRun")
     })
   })
   
   #updates on UI elements ----------
   
   updateReduction <- function()
   {
     session$sendCustomMessage("handler_disableAllButtons", "umapConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first and then re-run UMAP or tSNE or Diffusion Map above.")
       else {
         showModal(modalDialog(div('Analysis in Progress. This operation may take several minutes, please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("umapPlot_loader")
         
         #get input
         dims <- as.numeric(input$umapDimensions)
         type <- input$umapType
         
         #prepare metadata
         meta <- seurat_object@meta.data
         meta$Cell_id <- rownames(meta)
         meta <- meta[, ]
         reduc_data <- data.frame()
         
         #prepare colors
         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy])))
         
         #for all reductions
         seurat_object_reduc <- as.data.frame(seurat_object@reductions[[input$umapType]]@cell.embeddings)
         seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
         seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
         reduc_data <- left_join(seurat_object_reduc, meta)
         
         if(type == "umap" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="umap_1", y="umap_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder))+
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="UMAP 1", y="UMAP 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})
         }
         else if(type == "umap" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~umap_1, y=~umap_2, z=~UMAP_3, type="scatter3d", alpha = as.numeric(input$umapDotOpacity), mode="markers", color=as.formula(paste0('~', input$umapColorBy)),
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) )
           
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "tsne" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="tSNE_1", y="tSNE_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="tSNE 1", y="tSNE 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "tsne" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~tSNE_1, y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "dfm" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="DC_1", y="DC_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="DC 1", y="DC 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "dfm" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~DC_1, y=~DC_2, z=~DC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "pca" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="PC_1", y="PC_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="PC 1", y="PC 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "pca" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~PC_1, y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "phate" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="PHATE_1", y="PHATE_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="PHATE 1", y="PHATE 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "phate" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~PHATE_1, y=~PHATE_2, z=~PHATE_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else
         {
           session$sendCustomMessage("handler_alert", "There was an error with the number of dimensions of the plot. Please set Dimensions to 2D and update the plot.")
         }
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with drawing the resutls.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapConfirm")
     })
   }
   
   updateReduction2 <- function(typeR)
   {
     session$sendCustomMessage("handler_disableAllButtons", "umapConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first and then re-run UMAP or tSNE or Diffusion Map above.")
       else {
         showModal(modalDialog(div('Dimensional reduction in progress. This operation may take several minutes, please wait...')))
         shinyjs::show("umapPlot_loader")
         
         #get input
         dims <- as.numeric(input$umapDimensions)
         type <- typeR
         
         #prepare metadata
         meta <- seurat_object@meta.data
         meta$Cell_id <- rownames(meta)
         meta <- meta[, ]
         reduc_data <- data.frame()
         
         #prepare colors
         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy])))
         
         #for all reductions
         seurat_object_reduc <- as.data.frame(seurat_object@reductions[[type]]@cell.embeddings)
         seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
         seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
         reduc_data <- left_join(seurat_object_reduc, meta)
         print(head(reduc_data)) 
         
         if(type == "umap" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="umap_1", y="umap_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder))+
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="UMAP 1", y="UMAP 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})
         }
         else if(type == "umap" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~umap_1, y=~umap_2, z=~UMAP_3, type="scatter3d", alpha = as.numeric(input$umapDotOpacity), mode="markers", color=as.formula(paste0('~', input$umapColorBy)),
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) )
           
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "tsne" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="tSNE_1", y="tSNE_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="tSNE 1", y="tSNE 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "tsne" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~tSNE_1, y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "dfm" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="DC_1", y="DC_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="DC 1", y="DC 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "dfm" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~DC_1, y=~DC_2, z=~DC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "pca" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="PC_1", y="PC_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="PC 1", y="PC 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "pca" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~PC_1, y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else if(type == "phate" & dims == 2 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- ggplot(data=reduc_data, aes_string(x="PHATE_1", y="PHATE_2", fill=input$umapColorBy)) +
             geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
             scale_fill_manual(values = cols)+
             scale_size()+
             theme_bw() +
             theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                   axis.title.y = element_text(face = "bold", color = "black", size = 25),
                   axis.title.x = element_text(face = "bold", color = "black", size = 25),
                   legend.text = element_text(face = "bold", color = "black", size = 9),
                   legend.title = element_text(face = "bold", color = "black", size = 9),
                   legend.position="right",
                   title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
             labs(x="PHATE 1", y="PHATE 2", color="Cell type", title = "", fill="Color")
           output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
         }
         else if(type == "phate" & dims == 3 & dims <= ncol(seurat_object_reduc)-1)
         {
           p <- plot_ly(reduc_data, x=~PHATE_1, y=~PHATE_2, z=~PHATE_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                        marker = list(size = as.numeric(input$umapDotSize), 
                                      line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                        ),
                        colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
           output$umapPlot <- renderPlotly({print(p)})
         }
         else
         {
           session$sendCustomMessage("handler_alert", "There was an error with the number of dimensions of the plot. Please set Dimensions to 2D and update the plot.")
         }
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with drawing the resutls.")
     }, finally = {
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "umapConfirm")
     })
   }
   
   updateFeaturePair <- function()
   {
     if (!identical(seurat_object, NULL) & input$findMarkersFeaturePairReductionType != "-" & input$findMarkersFeaturePair1 %in% rownames(seurat_object) & input$findMarkersFeaturePair2 %in% rownames(seurat_object))
     {
       geneS1 <- input$findMarkersFeaturePair1
       geneS2 <- input$findMarkersFeaturePair2
       label_x <- ""
       label_y <- ""
       show_label <- as.logical(input$findMarkersFeaturePairLabels)
       order_exp <- as.logical(input$findMarkersFeaturePairOrder)
       minq <- paste0("q", input$findMarkersFeaturePairMinCutoff)
       maxq <- paste0("q", input$findMarkersFeaturePairMaxCutoff)
       blendThr <- as.numeric(input$findMarkersBlendThreshold)
       
       if(input$findMarkersFeaturePairReductionType == "umap")
       {
         label_x <- "umap_1"
         label_y <- "umap_2"
       }
       else if(input$findMarkersFeaturePairReductionType == "tsne")
       {
         label_x <- "tSNE_1"
         label_y <- "tSNE_2"
       }
       else if(input$findMarkersFeaturePairReductionType == "dfm")
       {
         label_x <- "DC_1"
         label_y <- "DC_2"
       }
       else if(input$findMarkersFeaturePairReductionType == "pca")
       {
         label_x <- "PC_1"
         label_y <- "PC_2"
       }
       else if(input$findMarkersFeaturePairReductionType == "phate")
       {
         label_x <- "PHATE_1"
         label_y <- "PHATE_2"
       }
       
       plot_temp <- FeaturePlot(seurat_object, features = c(geneS1, geneS2), blend.threshold = blendThr, 
                                pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red", "dodgerblue4"), 
                                order = order_exp, reduction = input$findMarkersFeaturePairReductionType, blend = TRUE, max.cutoff = maxq, min.cutoff = minq)
       
       plot <- FeaturePlot(seurat_object, features = c(geneS1, geneS2), blend.threshold = blendThr, 
                           pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red", "dodgerblue4"), raster = F,
                           order = order_exp, reduction = input$findMarkersFeaturePairReductionType, blend = TRUE, max.cutoff = maxq, min.cutoff = minq) +
         theme_bw() +
         theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
               axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
               axis.title.y = element_text(face = "bold", color = "black", size = 25),
               axis.title.x = element_text(face = "bold", color = "black", size = 25),
               legend.text = element_text(face = "bold", color = "black", size = 9),
               legend.title = element_text(face = "bold", color = "black", size = 9),
               legend.position="none",
               title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
         labs(x=label_x, y=label_y, color="Normalized\nexpression")
       gp1 <- plotly::ggplotly(plot[[1]] + xlim(min(plot_temp[[1]]$data[label_x]), max(plot_temp[[1]]$data[label_x])) + ylim(min(plot_temp[[1]]$data[label_y]), max(plot_temp[[1]]$data[label_y])) )
       gp2 <- plotly::ggplotly(plot[[2]] + xlim(min(plot_temp[[1]]$data[label_x]), max(plot_temp[[1]]$data[label_x])) + ylim(min(plot_temp[[1]]$data[label_y]), max(plot_temp[[1]]$data[label_y])) )
       gp3 <- plotly::ggplotly(plot[[3]] + xlim(min(plot_temp[[1]]$data[label_x]), max(plot_temp[[1]]$data[label_x])) + ylim(min(plot_temp[[1]]$data[label_y]), max(plot_temp[[1]]$data[label_y])) )
       gp4 <- plotly::ggplotly(plot[[4]]+theme(title = element_text(face = "bold", color = "black", size = 15)))
       
       output$findMarkersFPfeature1 <- renderPlotly({ gp1 })
       output$findMarkersFPfeature2 <- renderPlotly({ gp2 })
       output$findMarkersFPfeature1_2 <- renderPlotly({ gp3 })
       output$findMarkersFPcolorbox <- renderPlotly({ gp4 })  
     }
   }
   
   updateClusterTab <- function()
   {
     cluster_df <- as.data.frame(table(seurat_object$seurat_clusters)) #as.data.frame(table(Idents(seurat_object)))
     colnames(cluster_df)[1] <- "Cluster"
     colnames(cluster_df)[2] <- "Number of cells"
     cluster_df$`% of cells per cluster` <-  as.double(sprintf("%.3f", as.double( cluster_df$`Number of cells`/length(seurat_object@meta.data$orig.ident)*100 )))
     output$clusterTable <- renderDataTable(cluster_df, options = list(pageLength = 10))
     export_clustertable_RNA <<- cluster_df
   }
  
   updateUtilitiesAssays <- function()
   {
     getAssays <- names(seurat_object@assays)
     updateSelectInput(session, "utilitiesActiveAssay", choices = getAssays, selected = DefaultAssay(seurat_object))
   }
   # 
   updateQC_choices <- function()
   {
     updateSliderInput(session, "minUniqueGenes", min = min(init_seurat_object$nFeature_RNA)-1, max = max(init_seurat_object$nFeature_RNA)-2, value = as.numeric(quantile(init_seurat_object$nFeature_RNA, probs = 0.02)))
     updateSliderInput(session, "maxUniqueGenes", min = min(init_seurat_object$nFeature_RNA)+1, max = max(init_seurat_object$nFeature_RNA)+1, value = as.numeric(quantile(init_seurat_object$nFeature_RNA, probs = 0.98)))
   }
   # 
   updateMaxHVGs <- function()
   {
     if(length(rownames(seurat_object)) < 8000)
       updateSliderInput(session = session, inputId = "nHVGs", min = 200, max = length(rownames(seurat_object)), value = min(2000, length(rownames(seurat_object))))
   }
   # 
   updatePCs_selection <- function(numPCs)
   {
     updateSliderInput(session = session, inputId = "snnPCs", value = numPCs)
     updateSliderInput(session = session, inputId = "umapPCs", value = numPCs)
     updateSliderInput(session = session, inputId = "doubletsPCs", value = numPCs)
   }
   # 
   updateUmapTypeChoices <- function(type)
   {
     reductions_choices <<- c(reductions_choices, type)
     updateSelectInput(session, "umapType", choices = reductions_choices, selected = type)
     updateSelectInput(session, "findMarkersReductionType", choices = reductions_choices, selected = type)
     updateSelectInput(session, "findMarkersFeaturePairReductionType", choices = reductions_choices, selected = type)
     updateSelectInput(session, "doubletsReduction", choices = reductions_choices, selected = type)
     updateSelectInput(session, "cellCycleReduction", choices = reductions_choices, selected = type)
     updateSelectInput(session, "trajectoryReduction", choices = reductions_choices, selected = type)
   }
   # #update metadata RNA and export_RNA_table
   updateMetadata <- function()
   {
     seurat_object@meta.data$Cell_id <<- rownames(seurat_object@meta.data)
     #seurat_object@meta.data$percent.mt <<- as.double(sprintf("%.3f", as.double(seurat_object@meta.data$percent.mt)))

     #--3 decimal points
     for(i in 1:length(colnames(seurat_object@meta.data)))
     {
       if(is.double(seurat_object@meta.data[, i]))
       {
         seurat_object@meta.data[, i] <- as.double(sprintf("%.3f", as.double(seurat_object@meta.data[, i])))
       }
     }
     #--

     output$metadataTable <- renderDataTable(seurat_object@meta.data, options = list(pageLength = 10))
     export_metadata_RNA <<- seurat_object@meta.data
     updateRegressOut()
   }
   # #function update selectInput
   updateSelInpColor <- function()
   {
     factors <-colnames(seurat_object@meta.data)
     # tableMeta <- seurat_object@meta.data
     # f <- sapply(tableMeta, is.factor)
     # factors <- c(colnames(tableMeta[, f]),"customclassif")
     updateSelectInput(session, "utilitiesActiveClusters", choices = factors)
     updateSelectInput(session, "umapColorBy", choices = c(factors,"customclassif"))
     updateSelectInput(session, "samplecolumn", choices = factors)
     updateSelectInput(session, "qcColorBy", choices = factors)
     updateSelectInput(session, "utilitiesActiveAssay", choices = factors)
     updateSelectInput(session, "clusterGroupBy", choices = factors)
     updateSelectInput(session, "annotateClustersConfirmtype", choices = c(factors,"customclassif"))
   }
   updateRegressOut <- function()
   {
     variables <- colnames(seurat_object@meta.data)
     updateSelectInput(session, "normalizeRegressColumns", choices = variables)
   }

   updateGeneSearchFP <- function()
   {
     if(DefaultAssay(seurat_object) == "SCT")
       total_genes <- rownames(seurat_object@assays$SCT$counts)
     else if(DefaultAssay(seurat_object) == "tfactivity")
       total_genes <- rownames(seurat_object@assays$tfactivity$counts)
     else if(DefaultAssay(seurat_object) == "RNA")
       total_genes <- rownames(seurat_object@assays$RNA$counts)
     else if(DefaultAssay(seurat_object) == "integrated")
       total_genes <- rownames(seurat_object@assays$integrated$counts)

     #get signature names or numeric columns
     colnames(seurat_object@meta.data)
     tableMeta <- seurat_object@meta.data
     f <- sapply(tableMeta, is.numeric)
     sig_names <- colnames(tableMeta[, f])

     updateSelectizeInput(session, 'findMarkersGeneSelect', choices = c(total_genes, sig_names), server = TRUE) # server-side selectize drastically improves performance
     updateSelectizeInput(session, 'findMarkersGeneSelect2', choices = c(total_genes, sig_names), server = TRUE)
     updateSelectizeInput(session, 'findMarkersFeaturePair1', choices = total_genes, server = TRUE)
     updateSelectizeInput(session, 'findMarkersFeaturePair2', choices = total_genes, server = TRUE)
   }
   
   
   # This function is called after a new input file has been uploaded
   # and is responsible for clearing all generated plots across all tabs
   # @param fromDataInput: If TRUE, clears all, including QC, else skips clearing QC and metadata table
   cleanAllPlots <- function(fromDataInput){
     # renderPlotly
     if (fromDataInput){
       output$nFeatureViolin <- NULL
       output$totalCountsViolin <- NULL
       output$mitoViolin <- NULL
       output$mtCounts <- NULL
       output$genesCounts <- NULL
       output$cellStats <- NULL
       output$filteredNFeatureViolin <- NULL
       output$filteredTotalCountsViolin <- NULL
       output$filteredMitoViolin <- NULL
       output$filteredMtCounts <- NULL
       output$filteredGenesCounts <- NULL
       output$filteredCellStats <- NULL
       output$metadataTable <- NULL
       output$grnMatrixRNA <- NULL
       user_dir_pyscenic <<- ""
       user_dir <<- ""
     }
     output$elbowPlotPCA <- NULL
     output$PCAscatter <- NULL
     output$PCAloadings <- NULL
     output$PCAheatmap <- NULL
     output$gProfilerManhattan <- NULL
     output$annotateClustersCIPRDotplot <- NULL
     output$annotateClustersTissueplot <-NULL
     output$annotateClustersBubbleplot<- NULL
     output$ligandReceptorFullHeatmap <- NULL
     output$ligandReceptorCuratedHeatmap <- NULL
     output$hvgScatter <- NULL
     output$hvgTop10Stats <- NULL
     output$cellCyclePCA <- NULL
     output$cellCycleBarplot <- NULL
     output$clusterBarplot <- NULL
     output$umapPlot <- NULL
     output$findMarkersHeatmap <- NULL
     output$findMarkersDotplot <- NULL
     output$findMarkersFeaturePlot <- NULL
     output$findMarkersFPcolorbox <- NULL
     output$findMarkersFPfeature1_2 <- NULL
     output$findMarkersFPfeature1 <- NULL
     output$findMarkersFPfeature2 <- NULL
     output$findMarkersViolinPlot <- NULL
     output$findMarkersVolcanoPlot <- NULL
     output$grnHeatmapRNA <- NULL
     output$grnHeatmapRNA_DecoupleR <- NULL
     
     # renderPlot
     output$trajectoryPlot <- NULL
     output$trajectoryPseudotimePlot <- NULL
     
     # renderDataTable
     output$clusterTable <- NULL
     output$gProfilerTable <- NULL
     output$annotateClustersCIPRTable <- NULL
     output$findMarkersTable <- NULL
     output$grnMatrixRNA_DecoupleR <- NULL
     output$grnMatrixRNA <- NULL
     
     # renderPrint
     if (fromDataInput) output$cellStats <- NULL
     output$trajectoryText <- NULL
     
     # dittoDimPlot
     plot_D <- NULL
     plot_P <- NULL
     
     reductions_choices <<- c("-")
     updateSelectInput(session, "umapType", choices = reductions_choices)
     updateSelectInput(session, "findMarkersReductionType", choices = reductions_choices)
     updateSelectInput(session, "findMarkersFeaturePairReductionType", choices = reductions_choices)
     updateSelectInput(session, "doubletsReduction", choices = reductions_choices)
     updateSelectInput(session, "cellCycleReduction", choices = reductions_choices)
     updateSelectInput(session, "trajectoryReduction", choices = reductions_choices)
     
     #reset export table values
     export_metadata_RNA <<- ""
     export_loadingScoresTable_RNA <<- ""
     export_clustertable_RNA <<- ""
     export_markerGenes_RNA <<- ""
     export_enrichedTerms_RNA <<- ""
     export_annotation_RNA <<- ""
     export_ligandReceptor_full_RNA <<- ""
     export_ligandReceptor_short_RNA <<- ""
     export_scenicAUC_full_RNA <<- ""
     export_decoupleRZscores_full_RNA <- ""
   }
   # 
   updateInputLineageList <- function(lin_names)
   {
     updateSelectInput(session, "trajectoryLineageSelect", choices = lin_names)
   }

   updateInputLRclusters <- function() #update after clustering clusternames used in --> Volcano, gProfiler, slingshot, nichenetR
   {
     all_cluster_names <- (levels(seurat_object@meta.data[, 'seurat_clusters']))
     updateSelectInput(session, "ligandReceptorSender", choices = all_cluster_names)
     updateSelectInput(session, "ligandReceptorReciever", choices = all_cluster_names)
     updateSelectInput(session, "trajectoryStart", choices = all_cluster_names)
     updateSelectInput(session, "trajectoryEnd", choices = all_cluster_names)
     updateSelectInput(session, "findMarkersClusterSelect", choices = all_cluster_names)
     updateSelectInput(session, "gProfilerList", choices = all_cluster_names)
     updateSelectizeInput(session, "gProfilerFlameSelection", choices = all_cluster_names)
     updateSelectInput(session, "utilitiesDeleteCluster", choices = all_cluster_names)
     updateSelectInput(session, "utilitiesRenameOldName", choices = all_cluster_names)
   }
   
   # 
   updateFilteredQCplots <- function()
   {
     output$filteredNFeatureViolin <- renderPlotly(
       {
         p <- VlnPlot(seurat_object, features = c("nFeature_RNA"), pt.size = 0, group.by = "orig.ident",
                      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'orig.ident'])))) +
           theme_bw() +
           geom_hline(yintercept=c(as.numeric(input$minUniqueGenes), as.numeric(input$maxUniqueGenes)), linetype="dashed", color = "red", size=1) +
           theme(
             plot.title = element_blank(),
             axis.title.x = element_blank(),
             legend.position = "none") +
           labs(title = "", y="Genes detected/cell")
         plotly::ggplotly(p, tooltip = c("x", "y"))
       }
     )

     output$filteredTotalCountsViolin <- renderPlotly(
       {
         p <- VlnPlot(seurat_object, features = c("nCount_RNA"), pt.size = 0, group.by = "orig.ident",
                      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'orig.ident'])))) +
           theme_bw() +
           theme(
             plot.title = element_blank(),
             axis.title.x = element_blank(),
             legend.position = "none")+
           labs(title = "", y="Total counts/cell")
         plotly::ggplotly(p, tooltip = c("x", "y"))
       }
     )
   #   
     output$filteredMitoViolin <- renderPlotly(
       {
         p <- VlnPlot(seurat_object, features = c("percent.mt"), pt.size = 0, group.by = "orig.ident",
                      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'orig.ident'])))) +
           theme_bw() +
           geom_hline(yintercept= as.numeric(input$maxMtReads), linetype="dashed", color = "red", size=1) +
           theme(
             plot.title = element_blank(),
             axis.title.x = element_blank(),
             legend.position = "none")+
           labs(title = "", y="% of reads mapped to mitochondrial genome/cell")
         plotly::ggplotly(p, tooltip = c("x", "y"))
       }
     )

     output$filteredMtCounts <- renderPlotly(
       {
         p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", raster = F,
                             cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'orig.ident']))))
         gp <- plotly::ggplotly(p)
         print(gp)
       }
     )

     output$filteredGenesCounts <- renderPlotly(
       {
         p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", raster = F,
                             cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'orig.ident']))))
         gp <- plotly::ggplotly(p)
         print(gp)
       }
     )

     output$filteredCellStats <- renderPrint(
       {
         cat(paste0("Total number of cells: ", nrow(seurat_object@meta.data)))
       }
     )
   }

   # # Helper Functions ####

   sendToAPI <- function(apiLink, jsonBody) {
     response <- httr::POST(apiLink, body = jsonBody, encode = "json")
     urlFromResponse <- httr::content(response, as = "parsed")$url
     session$sendCustomMessage("handler_browseUrl", urlFromResponse)
   }

 
   # 
   # ###disable Tabs
   # 
   # disableTabsRNA <- function()
   # {
   #   hideTab(inputId="uploadTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="qcTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="pcaTabPanel", target="scRNA-seq: PCA", session = session)
   #   hideTab(inputId="clusteringTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="umapTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="findMarkersTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="featuresTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId ="doubletDetectionTabPanel", target = "scRNA-seq", session = session)
   #   hideTab(inputId="gProfilerTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="annotateClustersTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="trajectoryTabPanel", target="scRNA-seq", session = session)
   #   hideTab(inputId="grnTabPanel", target="scRNA-seq", session = session)
   #   
   #   
   #   ############################ show ATAC and make it selected
   #   
   #   showTab(inputId="uploadTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="qcTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="pcaTabPanel", target="scATAC-seq: LSI", session = session)
   #   showTab(inputId="clusteringTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="umapTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="findMarkersTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="featuresTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="doubletDetectionTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="gProfilerTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="annotateClustersTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="trajectoryTabPanel", target="scATAC-seq", session = session)
   #   showTab(inputId="grnTabPanel", target="scATAC-seq", session = session)
   #   
   #   updateTabsetPanel(inputId="uploadTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="qcTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="pcaTabPanel", selected="scATAC-seq: LSI", session = session)
   #   updateTabsetPanel(inputId="clusteringTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="umapTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="findMarkersTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="featuresTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="doubletDetectionTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId = "ATAC_markers_tabs", selected = "Marker genes (ATAC)", session = session)
   #   updateTabsetPanel(inputId="gProfilerTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="annotateClustersTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="trajectoryTabPanel", selected="scATAC-seq", session = session)
   #   updateTabsetPanel(inputId="grnTabPanel", selected="scATAC-seq", session = session)
   # }
   # 
   # disableTabsATAC <- function()
   # {
   #   hideTab(inputId="uploadTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="qcTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="pcaTabPanel", target="scATAC-seq: LSI", session = session)
   #   hideTab(inputId="clusteringTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="umapTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="findMarkersTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="featuresTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId ="doubletDetectionTabPanel", target = "scATAC-seq", session = session)
   #   hideTab(inputId="gProfilerTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="annotateClustersTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="trajectoryTabPanel", target="scATAC-seq", session = session)
   #   hideTab(inputId="grnTabPanel", target="scATAC-seq", session = session)
   # 
   #   ############################ show RNA and make it selected
   # 
   #   showTab(inputId="uploadTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="qcTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="pcaTabPanel", target="scRNA-seq: PCA", session = session)
   #   showTab(inputId="clusteringTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="umapTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="findMarkersTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="featuresTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId = "doubletDetectionTabPanel", target = "scRNA-seq", session = session)
   #   showTab(inputId="gProfilerTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="annotateClustersTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="trajectoryTabPanel", target="scRNA-seq", session = session)
   #   showTab(inputId="grnTabPanel", target="scRNA-seq", session = session)
   # 
   #   updateTabsetPanel(inputId="uploadTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="qcTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(session, inputId = "qc_tabs_rna", selected = "Pre-filtering plots")
   #   updateTabsetPanel(inputId="pcaTabPanel", selected="scRNA-seq: PCA", session = session)
   #   updateTabsetPanel(inputId="clusteringTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="umapTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="findMarkersTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="featuresTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="doubletDetectionTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="gProfilerTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="annotateClustersTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="trajectoryTabPanel", selected="scRNA-seq", session = session)
   #   updateTabsetPanel(inputId="grnTabPanel", selected="scRNA-seq", session = session)
   # }
   
   observeEvent(input$annotateClustersConfirmtype, {
     start.time <- Sys.time()
     print(paste0("!!!!-----scType_start-----", start.time))
     session$sendCustomMessage("handler_disableAllButtons", "annotateClustersConfirmtype")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
          else {
          showModal(modalDialog(div('Cluster annotation in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("annotateClustersCIPRDotplot_loader")
         es.max = sctype_score(scRNAseqData = seurat_object[["SCT"]]@scale.data, scaled = TRUE, 
                               gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
         cL_results = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
           es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
           head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
         }))
         sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
         sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
         for(j in unique(sctype_scores$cluster)){
           cl_type = sctype_scores[sctype_scores$cluster==j,]; 
           seurat_object@meta.data$customclassif[seurat_object@meta.data$seurat_clusters == j] <<- as.character(cl_type$type[1])
         }
         print(head(seurat_object@meta.data))
         updateMetadata()
         updateSelInpColor()
         tissue_guess <-auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seurat_object, scaled = TRUE, assay = "SCT")
         p<-ggplot(tissue_guess, aes(x=reorder(tissue,-score),y=score,fill=rgb(0.8,0.1,0.1,0.6))) +
           theme_bw()+
           geom_bar(stat="identity")+
           theme(legend.position = "none",axis.text.x = element_text(vjust=0.5, hjust=1, angle = 90))+
           scale_size(range = c(7, 7), guide = "none") +
           xlab("Tissue") +ylab("Summary score")
         output$annotateClustersTissueplot <- renderPlotly({ print(p)})
         cL_results<-cL_results[!duplicated(cL_results),]
         cL_results1=cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
                  # prepare nodes
         nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
         for (i in 1:length(unique(cL_results1$cluster))){
           dt_tmp = cL_results1[cL_results1$cluster == unique(cL_results1$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
         }
         nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
         files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
         nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
         nodes<-subset(nodes,shortName!="ECs")
         mygraph <- graph_from_data_frame(edges, vertices=nodes)

         # Make the graph
         gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
           geom_node_circle(aes(filter=ord==1), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour)), alpha=0.9) +
           theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, size = I(3), parse = T), segment.linetype="dotted")
         output$annotateClustersBubbleplot <- renderPlotly({ print(gggr)})
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with cluster annotation.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----scType_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "annotateClustersConfirm")
     })
     })
   observeEvent(input$annotateClustersConfirm, {
     start.time <- Sys.time()
     print(paste0("!!!!-----CIPR_start-----", start.time))
     
     session$sendCustomMessage("handler_disableAllButtons", "annotateClustersConfirm")
     tryCatch({
       if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data first.")
       else if (identical(seurat_object@misc$markers, NULL)) session$sendCustomMessage("handler_alert", "Please, execute the gene differential analysis at the MARKERS' IDENTIFICATION tab first.")
       else {
         showModal(modalDialog(div('Cluster annotation in progress. Please wait...'))) #position:absolute;top:50%;left:50%
         shinyjs::show("annotateClustersCIPRDotplot_loader")
         
         marker_genes <- seurat_object@misc$markers
         
         avgexp <- AverageExpression(seurat_object)
         avgexp <- as.data.frame(avgexp$SCT)
         avgexp$gene <- rownames(avgexp)
         input_CIPR <- marker_genes
         
         #Average expresssion methods
         if(input$annotateClustersMethod %in% c("all_genes_spearman", "all_genes_pearson"))
         {
           input_CIPR <- avgexp
         }
         
         CIPR(input_dat = input_CIPR,
              comp_method = input$annotateClustersMethod, 
              reference = input$annotateClustersReference, 
              keep_top_var = as.numeric(input$annotateClustersSlider),
              plot_ind = F,
              plot_top = F)
         
         CIPR_top_results$index <- as.character(CIPR_top_results$index)
         CIPR_top_results$index <- factor(CIPR_top_results$index, levels = as.character(seq(1:length(seurat_object$seurat_clusters))))
         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(CIPR_top_results$cluster)))
         
         p <- ggplot(CIPR_top_results, aes(x=index, y=identity_score, fill = cluster, size=7)) +
           theme_bw() +
           geom_point(alpha=1, shape=21) +
           scale_size(range = c(7, 7), guide = "none") + 
           scale_x_discrete(labels=CIPR_top_results$long_name) +
           scale_fill_manual(values = cols) +
           labs(x="") + 
           theme(legend.position="top",
                 axis.text.x = element_text(vjust=0.5, hjust=1, angle = 90))
         
         output$annotateClustersCIPRTable <- renderDataTable(CIPR_top_results[], options = list(pageLength = 20), rownames = F) #remove Description
         export_annotation_RNA <<- CIPR_top_results[]
         
         output$annotateClustersCIPRDotplot <- renderPlotly({ print(p)})
       }
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "There was an error with cluster annotation.")
     }, finally = {
       end.time <- Sys.time()
       time.taken <- as.numeric (end.time - start.time, units = "mins")
       print(paste0("!!!!-----CIPR_end_Mins_taken-----", time.taken))
       
       removeModal()
       Sys.sleep(1)
       session$sendCustomMessage("handler_enableAllButtons", "annotateClustersConfirm")
     })
   })
   
   output$annotationRNAExport <- downloadHandler(
     filename = function() { 
       paste("clusterAnnotationTableRNA-", Sys.Date(), ".txt", sep="")
     },
     content = function(file) {
       write.table(export_annotation_RNA, file, sep = "\t", quote = F, row.names = F)
     })
   #------------------Utilities------------------------------------------
   #save rna object as rds file
   output$utilitiesConfirmExport1 <- downloadHandler(
     filename = function() { 
       paste("processed_seurat_object-", Sys.Date(), ".RDS", sep="")
     },
     content = function(file) {
       saveRDS(seurat_object, file)
     })
   
   output$utilitiesConfirmExport2 <- downloadHandler(
     filename = function() { 
       paste("processed_seurat_object-", Sys.Date(), ".RDS", sep="")
     },
     content = function(file) {
       saveRDS(seurat_object, file)
     })
   
   output$utilitiesConfirmExport <- downloadHandler(
     filename = function() { 
       paste("processed_seurat_object-", Sys.Date(), ".RDS", sep="")
     },
     content = function(file) {
       saveRDS(seurat_object, file)
     })
   
   #rename-merge-delete clusters
   observeEvent(input$utilitiesConfirmRename, {
     if(!identical(seurat_object, NULL) & input$utilitiesRenameOldName != "-" & input$utilitiesRenameNewName != "")
     {
       showModal(modal_confirm_utilities_rename) 
     }
   })
   
   observeEvent(input$utilitiesConfirmDelete, {
     if(!identical(seurat_object, NULL) & input$utilitiesDeleteCluster != "-")
     {
       showModal(modal_confirm_utilities_delete)
     }
   })
   
   
   #alert2 rename cluster
   observeEvent(input$confirmUtilitiesRename, {
     old_idents_for_change <- seurat_object$seurat_clusters
     seurat_object$seurat_clusters <<- gsub(pattern = paste0("^",input$utilitiesRenameOldName,"$"), replacement = input$utilitiesRenameNewName, x = old_idents_for_change)
     seurat_object$seurat_clusters <<- as.factor(seurat_object$seurat_clusters)
     Idents(seurat_object) <<- seurat_object$seurat_clusters
     updateInputLRclusters()
     updateMetadata()
     #showNotification("Successful operation.", type = "message", closeButton = T)
     shinyalert(title = "Cluster rename operation completed succesfully!", type = "success")
     removeModal()
   })
   
   observeEvent(input$cancelUtilitiesRename, {
     removeModal()
   })
   
   #change active clusters
   observeEvent(input$utilitiesConfirmChangeCluster, {
     tryCatch({
       seurat_object$seurat_clusters <<- seurat_object@meta.data[, input$utilitiesActiveClusters]
       Idents(seurat_object) <<- seurat_object$seurat_clusters
       updateInputLRclusters()
       updateMetadata()
       updateClusterTab()
       shinyalert(title = "Active clustering column changed succesfully", type = "success")
     }, error = function(e) {
       print(paste("Error :  ", e))
       shinyalert(title = "There was an error with changing the clustering column.", type = "error")
     }, finally = {
       
     })
   })
   
   #change active assay 
   observeEvent(input$utilitiesConfirmChangeAssay, {
     tryCatch({
       DefaultAssay(seurat_object) <<- input$utilitiesActiveAssay
       updateGeneSearchFP()
       shinyalert(title = "Default assay was changed", type = "success")
     }, error = function(e) {
       print(paste("Error :  ", e))
       shinyalert(title = "There was an error with changing the default assay.", type = "error")
     }, finally = {
       
     })
     #print("def assay")
     #print(DefaultAssay(seurat_object))
   })
   
   observeEvent(input$utilitiesPlotCommands, {
     if(!is.null(seurat_object))
     {
       if(length(seurat_object@commands) > 0)
       {
         output$history <- renderPrint(
           {
             for(i in 1:length(seurat_object@commands))
             {
               print(seurat_object@commands[[i]])
               print("-------------------------------------------------------------")
             }
           })
       }
       else
       {
         Sys.sleep(1)
         shinyalert("Error. No command history was found.", type = "error")
       }
     }
     else
     {
       Sys.sleep(1)
       shinyalert("Error. No project is loaded.", type = "error")
     }
   })
   
   output$uploadMetadataExportRNA <- downloadHandler(
     filename = function() { 
       paste("metadataTableRNA-", Sys.Date(), ".txt", sep="")
     },
     content = function(file) {
       write.table(export_metadata_RNA, file, sep = "\t", quote = F, row.names = F)
     })

   #alert
   observeEvent(input$removeProject, {
     #remove RNAseq project
     init_seurat_object <<- NULL
     seurat_object <<- NULL
     cleanAllPlots(TRUE)
     
   hideAllLoaders()
     Sys.sleep(1)
     shinyalert(title = "Project removed succesfully! You can now upload a new project.", type = "success")
     removeModal()
   })
   
   observeEvent(input$cancelProjectRemoval, {
     removeModal()
   })
   modal_confirm <- modalDialog(
     "This operation cannot be completed, because another project is already loaded. Do you want to discard it?",
     title = "New project",
     footer = tagList(
       actionButton(inputId="removeProject", label="Yes"),
       actionButton("cancelProjectRemoval", "No")
     )
   )

   modal_confirm_utilities_rename <- modalDialog(
     "By renaming a cluster, you will need to update existing plots that have already been produced with the previous clustering scheme.
    By changing the name of a cluster to an already existing cluster name will lead to merging of the two clusters.
    Do you want to proceed?",
     title = "Cluster rename",
     footer = tagList(
       actionButton(inputId="confirmUtilitiesRename", label="Yes"),
       actionButton("cancelUtilitiesRename", "No")
     )
   )

   #alert3 delete cluster
   observeEvent(input$confirmUtilitiesDelete, {
     ident_to_remove <- input$utilitiesDeleteCluster
     seurat_object <<- subset(seurat_object, idents = ident_to_remove, invert = TRUE)
     seurat_object$seurat_clusters <<- Idents(seurat_object)
     updateInputLRclusters()
     updateMetadata()
     #showNotification("Successful operation.", type = "message", closeButton = T)
     shinyalert(title = "Cluster deletion operation completed succesfully", type = "success")
     removeModal()
   })

   observeEvent(input$cancelUtilitiesDelete, {
     removeModal()
   })

   modal_confirm_utilities_delete <- modalDialog(
     "By deleting a cluster, the user is advised to repeat various steps of the analysis and/or update existing plots.
    In more detail, the scaled expression values calculated in the tab \"DATA NORMALIZATION & SCALING\" should be recalculated as the
    final number of cell has changed.
    Do you want to proceed?",
     title = "Cluster deletion",
     footer = tagList(
       actionButton(inputId="confirmUtilitiesDelete", label="Yes"),
       actionButton("cancelUtilitiesDelete", "No")
     )
   )
}
)