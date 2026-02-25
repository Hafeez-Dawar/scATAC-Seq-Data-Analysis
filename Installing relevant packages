# Update BiocManager to the latest version for R 4.4
install.packages("BiocManager")
BiocManager::install(version = "3.19")

# Then install Signac
BiocManager::install("Signac")
BiocManager::install("hdf5r")


# Install CRAN packages first
install.packages("RCurl")
install.packages("XML")

# Install GenomicAlignments
BiocManager::install("GenomicAlignments")

# Install remaining Bioconductor packages
BiocManager::install(
  c("BiocIO", "restfulr", "GenomicFeatures", "AnnotationFilter", 
    "rtracklayer", "ProtGenerics", "ensembldb"),
  ask = FALSE,    # Don't ask about updates
  update = FALSE, # Don't update packages
  force = FALSE   # Don't force reinstall unless needed
)

# Install biovizBase and EnsDb.Hsapiens.v86
BiocManager::install("biovizBase")
BiocManager::install("EnsDb.Hsapiens.v86", force = TRUE, ask = FALSE)

# 3. Install plotly and virdislite from scratch
install.packages("viridisLite")
install.packages("plotly")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(biovizBase)            
library(EnsDb.Hsapiens.v86)
library(viridisLite)
library(plotly)


# Download scATAC-Seq data from 10X Genomics
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5

wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv

wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz

wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi



#or you can download it into your local machine by using the download.file() function in R
# 1. Filtered peak-barcode matrix (HDF5)
download.file(
  url = "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5",
  destfile = "10k_pbmc_ATACv2/filtered_peak_bc_matrix.h5",
  mode = "wb"
)

# 2. Single-cell metadata CSV
download.file(
  url = "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  destfile = "10k_pbmc_ATACv2/singlecell.csv",
  mode = "wb"
)

# 3. Fragments file
download.file(
  url = "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  destfile = "10k_pbmc_ATACv2/fragments.tsv.gz",
  mode = "wb"
)
