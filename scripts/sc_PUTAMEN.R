library(Seurat)



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
    if (length(rownames(cellSplits[[i]])) < 500){
      n = length(rownames(cellSplits[[i]]))
    }
    else {
      n = 500
    }
    cell.list <- WhichCells(cellSplits[[i]], downsample = n)
    cellSplits_red[i] <- seurat[, cell.list]
  }  
  
  seurat = merge(cellSplits_red[[1]], y = cellSplits_red[2:length(unique(seurat$cellType))] )    
  return(seurat)
}



seurat = readRDS('sc_PUTAMEN/local.rds')

seurat<- seurat[, sample(colnames(seurat), size = 10000, replace=F)]

unique(seurat$supercluster_term)

unique(seurat$cell_type)

#Idents(seurat) = seurat$cell_type

seurat$cellType = seurat$supercluster_term

seurat = get_300_cells_per_celltype(seurat)

seurat$cellType = seurat$cell_type

seurat = remove_rare_celltypes(seurat)

sc<- FindVariableFeatures(seurat)#, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType

DimPlot(seurat, label = TRUE)   + NoLegend() +ggtitle('oary')

# rename myofibroblast to fibroblast

pbmc$sampleID = pbmc$donor_id
pbmc$cellID = colnames(pbmc)

saveRDS(pbmc, 'sc_PUTAMEN//E-MTAB-5214_tabula_putamen_seurat.rds')
count(seurat300$cellType)


DefaultAssay(object = seurat300)
