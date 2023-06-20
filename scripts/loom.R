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
files = list.files('Raw/Split/', pattern = '*.loom')
regions = str_remove(files, '_5000.loom')
print(regions)
#tissues = c("lung", "blood", "pancreas", "liver", 'coronaryartery', 'cortexofkidney', 'EBV-transformedlymphocyte', 'endocervix', 'ectocervix', 'uterus')

for (region in regions){
#tissues = c("lung", "blood", "pancreas", "liver", 'coronaryartery', 'cortexofkidney', 'EBV-transformedlymphocyte', 'endocervix', 'ectocervix', 'uterus')
    seurat = LoadLoom(loom_file)

    # get the CL ids for cell types
    seurat = get_300_cells_per_celltype(seurat)
    CL_ids = factor()
    cell_types = (unique(seurat$cellType))
    for (cell_type in cell_types){
      print(cell_type)
      CL_ids = append(CL_ids , get_semantic_tag(cell_type, 'CL'))
    }
    
    seurat$cell_type_ontology_term_id <- mapvalues(seurat$cell_type, from =  cell_types, to = CL_ids)

    # only get the cells that have CL ids
    seurat <- seurat[, grepl("CL", seurat$ell_type_ontology_term_id)]
    # Save the Seurat object to a file
    saveRDS(seurat, file = paste0("Raw/Split/", region, '_seurat.rds'))
}