
library('Seurat')
library(dplyr)
library(plyr)


source('scripts/R_functions.R')

seurat = readRDS('Raw/FallopianTube.rds')

seurat$cellType = seurat$cell_type

seurat300 = get_300_cells_per_celltype(seurat)

seurat300$sampleID = seurat300$donor_id
seurat300$cellID = colnames(seurat300)

saveRDS(seurat300, 'Processed/fallopiantube_seurat.rds')
