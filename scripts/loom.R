library(scater)
library(SeuratDisk)
library(reticulate)
library(Seurat)
library(loomR)
library(dplyr)
library(plyr)
library(ggplot2)
library(Matrix)

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
remove_rare_celltypes = function(seurat){
  row_count <- as.data.frame(x = table(seurat$cellType))
  print(row_count)
  non_rare = row_count[row_count$Freq >100 ,]
  
  #sc$cellType_freq = mapvalues(df_cellType$sc.cellType, from = row_count$Var1, to = row_count$Freq)
  
  seurat <- subset(seurat, subset = cellType %in% paste(non_rare$Var1))
  return(seurat)
  
}

seurat = LoadLoom('deconvolution/lin_BRAIN/Spinal cord_5000.loom')


seurat = get_300_cells_per_celltype(seurat)
seurat = remove_rare_celltypes(seurat)

unique(seurat$cellType)

Idents(seurat) = seurat$cellType



sc<- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)



all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType

DimPlot(pbmc, label = TRUE)  +ggtitle('Cerebellum_8000')

# rename myofibroblast to fibroblast
#pbmc$cellType = mapvalues(pbmc$cellType, from = 'myofibroblast cell', to = 'fibroblast')
pbmc$cellID= pbmc$sampleID
pbmc$sampleID = colnames(pbmc)

# Save the Seurat object to a file
saveRDS(pbmc, file = "deconvolution/lin_BRAIN/E-MTAB-5214_tabula_C1segmentofcervicalspinalcord_seurat.rds")


