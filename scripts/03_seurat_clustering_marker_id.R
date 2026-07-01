#This is the script for the analyses in R for independently replicating the snRNA analysis in Shamsi, F., Piper, M., Ho, LL. et al. 
#Vascular smooth muscle-derived Trpv1+ progenitors are a source of cold-induced thermogenic adipocytes. Nat Metab 3, 485–495 (2021). https://doi.org/10.1038/s42255-021-00373-z
#This script handles the seurat clustering, subclustering and analysis of markers to identify specialised TRPV1+ APCs after analysis of the data in python using scverse tools.



# PIPELINE STAGE: 3 of 3 (see repo README for full pipeline table)
#   01_qc_filtering_h5ad_export.R -> 02_scvi_integration_umap.ipynb -> THIS SCRIPT



# Load packages
library(reticulate)
use_condaenv("/home/abhi/miniconda3/envs/sc", required = TRUE) # Specifying the python to be used to avoid dependency issues.
library(Seurat)
library(sceasy)
library(dplyr)
library(ggplot2)

# Convert h5ad to seurat
H_BAT_integrated.seurat <- sceasy::convertFormat("H_BAT_nF_integrated.h5ad",
from = "anndata", to = "seurat")

# Verify the imported seurat object
H_BAT_integrated.seurat
dim(H_BAT_integrated.seurat)
head(H_BAT_integrated.seurat@meta.data)
Reductions(H_BAT_integrated.seurat)
Assays(H_BAT_integrated.seurat)

# Save the seurat object
saveRDS(H_BAT_integrated.seurat, file = "H_BAT_integrated.seurat.rds")

# Normalize expression
H_BAT_integrated.seurat <- NormalizeData(H_BAT_integrated.seurat)

# Find variable genes
H_BAT_integrated.seurat <- FindVariableFeatures(H_BAT_integrated.seurat)

# Scale data
H_BAT_integrated.seurat <- ScaleData(H_BAT_integrated.seurat)

# Plot the imported UMAP
DimPlot(H_BAT_integrated.seurat, reduction = "umap", label = FALSE)

# Build the nearest neighbour graph
H_BAT_integrated.seurat <- FindNeighbors(H_BAT_integrated.seurat, reduction = "scVI", dims = 1:20)

# Cluster cells
H_BAT_integrated.seurat <- FindClusters(H_BAT_integrated.seurat, resolution = 0.5)

# Plot clusters
umap1 <- DimPlot(H_BAT_integrated.seurat, reduction = "umap", label = TRUE)
 
ggsave(
filename = "UMAP1.png",
plot = umap1,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)

# Identify cluster markers
markers <- FindAllMarkers(H_BAT_integrated.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Save markers
write.csv(markers, "H_BAT_cluster_markers.csv",row.names = FALSE)

# Visualise markers
markerplot1 <- FeaturePlot(
H_BAT_integrated.seurat,
reduction = "umap",
features = c(
"PDGFRA",
"ADIPOQ",
"TRPV1"))

ggsave(
filename = "markerplot1.png",
plot = markerplot1,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)


# Violin plot
vlnplot1 <- VlnPlot(
H_BAT_integrated.seurat,
features = c("PDGFRA", "ADIPOQ"))

ggsave(
filename = "violinplot1.png",
plot = vlnplot1,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)

# Identified clusters 2 and 3 as adipocytes and APCs respectively
# Subset adipocytes (ADIPOQ) and adipocyte progenitor cells (PDGFRA).
H_BAT_subset <- subset(H_BAT_integrated.seurat, idents = c("2","3"))

# Find neighbours and recluster the subset
H_BAT_subset <- FindNeighbors(H_BAT_subset, reduction = "scVI", dims = 1:20)

H_BAT_subset <- FindClusters(H_BAT_subset,resolution = 0.8)

# Plot the umap for the new clustering
umap2 <- DimPlot(H_BAT_subset, reduction = "umap", label = TRUE)

ggsave(
filename = "umap2.png",
plot = umap2,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)


markers_subset <- FindAllMarkers(H_BAT_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Markers for the subset
write.csv(markers_subset, "H_BAT_subset_markers.csv",row.names = FALSE)
markerplot2 <- FeaturePlot(
H_BAT_subset,
reduction = "umap",
features = c("PDGFRA", "ADIPOQ", "TRPV1"),
pt.size = 0.2)

ggsave(
filename = "markerplot2.png",
plot = markerplot2,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)


# Violin plot for the marker expression
vlnplot2 <- VlnPlot(
H_BAT_subset,
features = c("PDGFRA", "TRPV1", "ADIPOQ"))

ggsave(
filename = "violinplot2.png",
plot = vlnplot2,
path = "/home/abhi/sc/plots/",
width = 7,
height = 5,
dpi = 300)
