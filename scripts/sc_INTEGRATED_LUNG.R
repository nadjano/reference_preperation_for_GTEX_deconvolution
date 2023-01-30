library(Seurat)

library('Seurat')
library(dplyr)
library(plyr)
library(ggplot2)

remove_rare_celltypes = function(seurat){
  row_count <- as.data.frame(x = table(seurat$cellType))
  
  non_rare = row_count[row_count$Freq >100 ,]
  
  #sc$cellType_freq = mapvalues(df_cellType$sc.cellType, from = row_count$Var1, to = row_count$Freq)
  
  seurat <- subset(seurat, subset = cellType %in% paste(non_rare$Var1))
  return(seurat)
  
}
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



seurat = readRDS('sc_LUNG/integratedLungCellAtlas_30000.rds')



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

saveRDS(pbmc, 'sc_LUNG///E-MTAB-5214_tabula_lung_seurat.rds')

sc= readRDS('sc_LUNG///E-MTAB-5214_tabula_lung_seurat.rds')

unique(sc$cellType)
