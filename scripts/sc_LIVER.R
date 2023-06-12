library(Seurat)
library(ggplot2)
library("plyr")

source('scripts/R_functions.R')


seurat = readRDS('Raw/liver.rds')
seurat$cellType = seurat$BroadCellType
seurat = get_300_cells_per_celltype(seurat)
seurat = remove_rare_celltypes(seurat)


sc<- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType
#DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0("liver", ' (#: ', ncol(pbmc), ')' ))
#rename clusters

#sample some donors
Idents(pbmc) = pbmc$cellType
pbmc= RenameIdents(pbmc,  "T/NK-cell" = 'T.NK.cell',  "B-cell" = "B.cell")
pbmc$cellType = Idents(pbmc)

pbmc$sampleID = pbmc$donor_id
pbmc$cellID = colnames(pbmc)
saveRDS(pbmc, 'Processed/liver_seurat.rds')


