
library('Seurat')
library(dplyr)
library(plyr)


source('scripts/R_functions.R')


# seurat = readRDS('Raw/GutCellAtlas.rds')
# seurat  = subset(seurat, subset = tissue %in% c('small intestine', 'large intestine'))
seurat  = readRDS('Raw/GutCellAtlas_reduced.rds')
write.table(unique(seurat$cell_type),sep= '\t',file ='Raw/GutCellAtlas_cellTypes.tsv' )
#print(head(seurat@meta.data))
# unique(seurat$cell_type)

# Idents(seurat) = seurat$cell_type

# seurat$cellType = seurat$cell_type

# seurat300 = get_300_cells_per_celltype(seurat)
# seurat300 = remove_rare_celltypes(seurat300)

# seurat300$tissue = droplevels(x = seurat300$tissue)
# seurat300$cellType = droplevels(x = seurat300$cellType)

# unique(seurat300$cell_type)
# sc<- FindVariableFeatures(seurat300, selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(sc)
# pbmc <- ScaleData(sc, features = all.genes)

# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# pbmc <- RunUMAP(pbmc, dims = 1:10)

# Idents(pbmc) = pbmc$cellType



# # rename myofibroblast to fibroblast

# pbmc$cellType = mapvalues(pbmc$cellType, from = 'myofibroblast cell', to = 'fibroblast')
# pbmc$sampleID = pbmc$donor_id
# pbmc$cellID = colnames(pbmc)

# saveRDS(pbmc, 'Processed/sigmoidcolon_seurat.rds')
# saveRDS(pbmc, 'Processed/transversecolon_seurat.rds')
# saveRDS(pbmc, 'Processed/smallintestinePeyerspatch_seurat.rds')