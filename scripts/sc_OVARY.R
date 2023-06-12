
library('Seurat')
library(dplyr)
library(plyr)


source('scripts/R_functions.R')

seurat = readRDS('Raw/ovary.rds')

seurat@assays$RNA

unique(seurat$cell_type)

Idents(seurat) = seurat$cell_type

seurat$cellType = seurat$cell_type

seurat300 = get_300_cells_per_celltype(seurat)

seurat300 = remove_rare_celltypes(seurat300)

sc<- FindVariableFeatures(seurat300, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType


pbmc$sampleID = pbmc$donor_id
pbmc$cellID = colnames(pbmc)

saveRDS(pbmc, 'Processed/ovary_seurat.rds')
count(seurat300$cellType)
