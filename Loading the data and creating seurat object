#R Version (R version 4.4.3)
library(Signac) # version 1.16.0
library(Seurat) # 5.4.0
library(GenomicRanges) # 1.56.2
library(ggplot2) # 4.0.1
library(patchwork) # 1.3.2
library(hdf5r)
library(biovizBase)            
library(EnsDb.Hsapiens.v86)

list.files()

#load the .h5 file
counts <-Read10X_h5("10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
# check the dimension
dim(counts)
# Look into the first few rows and columns
counts[1:12, 1:4]



#load the .csv file
metadata <- read.csv("10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv", header = TRUE, row.names = 1)
# Check the first few rows and columns
metadata[1:6, 1:6]
View(metadata)



#create a chromatin assay
chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"),
fragments = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz", min.cells = 10, min.features = 200)


#create seurat object
pbmc <- CreateSeuratObject(counts = chrom_assay, meta.data = metadata, assay = "peaks")
pbmc
View(pbmc@meta.data)


#save the seurat object
saveRDS(pbmc, file = "pbmc.rds")

#Load the seurat object
pbmc <- readRDS("pbmc.rds")
