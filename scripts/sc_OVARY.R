
library('Seurat')
library(dplyr)
library(plyr)


source('scripts/R_functions.R')

# read in .rds downloaded from cellXgene
seurat = readRDS('Raw/ovary.rds')

seurat$cellType = seurat$cell_type

# downsample cells to max 300 per cell type
seurat300 = get_300_cells_per_celltype(seurat)

seurat300$sampleID = seurat300$donor_id
seurat300$cellID = colnames(seurat300)

saveRDS(seurat300, 'Processed/ovary_seurat.rds')

