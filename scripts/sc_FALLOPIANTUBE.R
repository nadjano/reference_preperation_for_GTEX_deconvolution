
library('Seurat')
library(dplyr)
library(plyr)
library(ggplot2)
get_300_cells_per_celltype <- function(seurat){
  
  #Split cells by type
  cellSplits <- SplitObject(seurat, split.by = "cellType")
  cellSplits_red = cellSplits
  
  for (i in 1:length(unique(seurat$cellType))){
    if (length(rownames(cellSplits[[i]])) < 300){
      n = length(rownames(cellSplits[[i]]))
    }
    else {
      n = 300
    }
    cell.list <- WhichCells(cellSplits[[i]], downsample = n)
    cellSplits_red[i] <- seurat[, cell.list]
  }  
  
  seurat = merge(cellSplits_red[[1]], y = cellSplits_red[2:length(unique(seurat$cellType))] )    
  return(seurat)
}



seurat = readRDS('sc_FALLOPIANTUBE/local.rds')

GetAssayData(seurat, slot = "counts")
GetAssayData(seurat, slot = "data")


unique(seurat$cell_type)

Idents(seurat) = seurat$cell_type

seurat$cellType = seurat$cell_type

seurat300 = get_300_cells_per_celltype(seurat)



sc<- FindVariableFeatures(seurat300, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType

DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle('fallopian tube')

# rename myofibroblast to fibroblast

pbmc$cellType = mapvalues(pbmc$cellType, from = 'myofibroblast cell', to = 'fibroblast')
pbmc$sampleID = pbmc$donor_id
pbmc$cellID = colnames(pbmc)

saveRDS(pbmc, 'sc_FALLOPIANTUBE/E-MTAB-5214_tabula_fallopiantube_seurat.rds')
