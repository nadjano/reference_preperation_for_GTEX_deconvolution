library(Seurat)

library('Seurat')
library(dplyr)
library(plyr)
library(ggplot2)

source('scripts/R_functions.R')


seurat = readRDS('Raw/IntegratedLungCellAtlas.rds')



unique(seurat$ann_level_3)

Idents(seurat) = seurat$ann_level_2

seurat$cellType = seurat$ann_level_2

seurat300 = get_300_cells_per_celltype(seurat)

seurat300 = remove_rare_celltypes(seurat300)

sc<- FindVariableFeatures(seurat300, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)



Idents(pbmc) = pbmc$cellType

DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle('lung')


# rename myofibroblast to fibroblast

pbmc$sampleID = pbmc$donor_id
pbmc$cellID = colnames(pbmc)

saveRDS(pbmc, 'Processed/lung_C0.rds')

