
library(Seurat)
library(ggplot2)
library(dplyr)

source('scripts/R_functions.R')

HCL = readRDS('Raw/HumanCellLandscapes.rds')
HCL$cellType = HCL$cell_type

tissues =c('stomach','thyroid gland' )

for (tissue in tissues){
  
  split_sc_per_tissue(HCL, tissue, gsub(" ", "", paste0('Processed/', tissue, '_seurat.rds')))

}


