
#This is the script for the analyses in R for independently replicating the snRNA analysis in Shamsi, F., Piper, M., Ho, LL. et al. 
#Vascular smooth muscle-derived Trpv1+ progenitors are a source of cold-induced thermogenic adipocytes. Nat Metab 3, 485–495 (2021). https://doi.org/10.1038/s42255-021-00373-z
#This script handles the initial quality check and filtering after which the data is processed in python using scverse tools.


# PIPELINE STAGE: 1 of 3 (see repo README for full pipeline table)
# [raw CellRanger data] -> THIS SCRIPT -> 02_scvi_integration_umap.ipynb -> 03_seurat_clustering_marker_id.R


# Load all required packages----
library(DropletUtils) # Used here to read in cellranger data and construct barcode ranks from the barcodeRanks function.
library(ggplot2) # Used here to plot the barcode rank plot.
library(sceasy) # Used here to convert the filtered data into a format used in python. 


# Setting up important objects in order to automate the process----
base_dir <- "/home/abhi/sc" # Base directory where all the sample folders to be read in are situated.
plot_dir <- file.path(base_dir, "Barcode_rank_plots1") # Directory where all the generated barcode rank plots will be saved.
sample_names <- c("H_BAT_F_1", "H_BAT_F_5", "H_BAT_F_6", "H_BAT_F_7", "H_BAT_F_8",
                  "H_BAT_nF_1", "H_BAT_nF_2", "H_BAT_nF_3", "H_BAT_nF_4") # Sample folder names

# Plotting function to produce barcode rank plots for each sample.----
plot_barcode_rank <- function(br, sample_name) {
  rank_df <- as.data.frame(br) # Converting the DropletUtils result into a data frame ggplot can use
  #rank_df <- rank_df[rank_df$total > 0, ] # Subsetting to remove barcodes with zero UMIs as they represent empty cells. 
  #log10(0) is infinity and this produces an error when plotted on logscale in ggplot.
  
  knee_value <- metadata(br)$knee 
  inflection_value <- metadata(br)$inflection
  
  
  
  min_rank <- min(rank_df$rank, na.rm = TRUE) # Parameter to anchor the two labels to the actual range of the data so that we dont encounter the problem of labels not appearing 
  max_rank <- max(rank_df$rank, na.rm = TRUE) # Parameter to anchor the two labels to the actual range of the data so that we dont encounter the problem of labels not appearing 
  
  ggplot(rank_df, aes(x = rank, y = total)) +
    geom_point(size = 0.5, alpha = 0.5, color = "#264653") +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(yintercept = knee_value, linetype = "dashed", color = "#2A6F97") +
    geom_hline(yintercept = inflection_value, linetype = "dashed", color = "#40916C") +
    annotate("text", x = max_rank, y = knee_value, label = "Knee",
             color = "#2A6F97", hjust = 1, vjust = -0.5, size = 3.5) +
    annotate("text", x = min_rank, y = inflection_value, label = "Inflection",
             color = "#40916C", hjust = 0, vjust = -0.5, size = 3.5) +
    labs(title = paste("Barcode rank plot_", sample_name),
         x = "Barcode rank (log10)",
         y = "Total UMI count (log10)") +
    theme_classic(base_size = 13)
}
# End of plotting function

sce_list <- list() # List to hold each sample's SingleCellExperiment object.
br_list  <- list() # List to hold each sample's barcodeRanks result.
 
# Loop to automate the qc of each samples----
for (sample_name in sample_names) {
  
  message("Processing ", sample_name, " ...")
  
  # 1. Read in this sample
  sample_path <- file.path(base_dir, sample_name)
  sce_list[[sample_name]] <- read10xCounts(sample_path)
  
  # 2. Calculate its barcode ranks
  br_list[[sample_name]] <- barcodeRanks(counts(sce_list[[sample_name]]))
  
  # 3. Create the barcode rank plot 
  p <- plot_barcode_rank(br_list[[sample_name]], sample_name)
  
  # 4. Save it into Barcode_rank_plots, named Barcode_rank_plot_<samplename>.png
  ggsave(
    filename = file.path(plot_dir, paste0("Barcode_rank_plot_", sample_name, ".png")),
    plot = p, width = 7, height = 5, dpi = 300
  )
}


passing_samples <- c("H_BAT_nF_1", "H_BAT_nF_2", "H_BAT_nF_3", "H_BAT_nF_4") # Samples selected after qc

# Filtering of data to keep only top 1000 genes based on UMI counts and remove doublets(cells with UMI>20000).----



# Keep only top 1000 genes based on UMI counts
filtered_list <- list() # List to keep only top 1000 genes based on UMI counts.

for (sample_name in passing_samples) {
  
  br  <- br_list[[sample_name]]
  sce <- sce_list[[sample_name]]
  
  ordered_barcodes <- order(br$total, decreasing = TRUE)
  top1000 <- ordered_barcodes[1:1000]
  
  filtered_list[[sample_name]] <- sce[, top1000]
} # For loop to filter only op 1000 genes based on UMI counts.




# Remove doublets(cells with UMI>20000)
clean_list <- list() # List to keep only single cells

for (sample_name in passing_samples) {
  
  filtered <- filtered_list[[sample_name]]
  umi_totals <- colSums(counts(filtered))
  
  clean_list[[sample_name]] <- filtered[, umi_totals <= 20000]
} # For loop to remove doublets.


#Conversion of data into anndata formata to be used in python
for (sample_name in passing_samples) {
  sceasy::convertFormat(
    clean_list[[sample_name]], from = "sce", to = "anndata",
    outFile = file.path(base_dir, paste0(sample_name, ".h5ad"))
  )
} # For loop for converting the filtered samples that passed qc to h5ad format to be used in python for next step of analysis.
