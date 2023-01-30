
library(Seurat)
library(ggplot2)
library(dplyr)

remove_rare_celltypes = function(seurat){
  row_count <- as.data.frame(x = table(seurat$cellType))
  print(row_count)
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
HCL = readRDS('deconvolution/HumanCellLandscapes/local.rds')
#stomach = subset(HCL, subset = tissue =='stomach')
#thyroid = subset(HCL, subset = tissue =='thyroid gland')
artery = subset(HCL, subset = tissue =='artery') 
seurat = artery
seurat$cellType = seurat$cell_type

seurat = get_300_cells_per_celltype(seurat)
seurat = remove_rare_celltypes(seurat)


sc<- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType
DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0("artery", ' (#: ', ncol(pbmc), ')' ))
#rename clusters
pbmc$sampleID = pbmc$donor_id
pbmc$CELLID = colnames(pbmc)
saveRDS(pbmc, 'deconvolution/HumanCellLandscapes/HCL_artery_seurat.rds')
#sample some donors
Idents(pbmc) = pbmc$cellType
pbmc= RenameIdents(pbmc,  "T/NK-cell" = 'T.NK.cell',  "B-cell" = "B.cell")
pbmc
