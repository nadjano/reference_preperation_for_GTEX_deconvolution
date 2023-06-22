#library(scater)
library(SeuratDisk)
#library(reticulate)
library(Seurat)
#library(loomR)
library(dplyr)
library(plyr)
library(ggplot2)
library(Matrix)
suppressMessages(library(stringr))
suppressMessages(library(httr))
suppressMessages(library(jsonlite))

source('scripts/R_functions.R')

# load loom file that was split into brain regions and downsamples to 
# 5000 cells per region
files = paste0('Split/', list.files('Split/', pattern = '*.loom'))

for (loom_file in files){

    seurat = LoadLoom(loom_file)
    seurat = seurat[, !is.na(seurat$cellType) ]
    region = str_remove(basename(loom_file), '_5000.loom')
    # get the CL ids for cell types
    print(region)
    seurat = get_300_cells_per_celltype(seurat)
    CL_ids = factor()
    cell_types = (unique(seurat$cellType))
    print(cell_types)
    for (cell_type in cell_types){
      print(cell_type)
      cell_type_name = gsub("\\.", " ", cell_type)
      CL_ids = append(CL_ids , get_semantic_tag(cell_type_name, '/CL'))
    }
    
    seurat$cell_type_ontology_term_id <- mapvalues(seurat$cellType, from =  cell_types, to = CL_ids)

    #only get the cells that have CL ids
    seurat <- seurat[, grepl("CL", seurat$cell_type_ontology_term_id)]
    # Save the Seurat object to a file
    saveRDS(seurat, file = paste0("Split/Linnarsson_2022_", region, '_seurat.rds'))
}