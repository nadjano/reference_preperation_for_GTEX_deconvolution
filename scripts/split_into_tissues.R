library(Seurat)
library(ggplot2)
library(dplyr)

# Run functions 
args <- commandArgs(trailingOnly = TRUE)
seurat_filename = args[1]

seurat = readRDS(seurat_filename)

print(head(seurat))

tissues = unique(seurat$tissue)

for (tissue in tissues){
    seurat_tissue <- seurat[, seurat$tissue == tissue]

    # Downsampling to get max 300 cells per celltype
    # this is done to reduce the size of the dataset
    cell_types <- unique(seurat_tissue$cell_type)
    seurat_downsampled = list()
    for (ct in cell_types) {
    cells <- colnames(subset(seurat_tissue, cell_type == ct))
    if (length(cells) > 300) {
        cells <- sample(cells, 300)
        print(cells)
    }
    # append the downsamples cell types
    seurat_downsampled[ct] <- seurat_tissue[,cells]
    }
    if (length(cell_types) > 1) {
        seurat_downsampled = merge(
                            seurat_downsampled[[1]],
                            y = seurat_downsampled[2:length(cell_types)]
                            )
    } else {
        seurat_downsampled = seurat_downsampled[[1]]
    }  
    
    # Save output
    saveRDS(seurat_downsampled, paste0('Raw/split/', gsub(' ', '_', tissue), 'seurat.rds'))
}
