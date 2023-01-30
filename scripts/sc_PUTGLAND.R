
library(Seurat)
library(ggplot2)
library("plyr")

seurat = Read10X('sc_PUT_GLAND/')
equal_celltype_per_sample <- function(seurat){
  
  donor = c('TSP1', 'TSP2', 'TSP3', 'TSP4', 'TSP5')
  seurat$sampleID = sample(donor, size = length(colnames(seurat)), replace=T, prob =c(0.2, 0.2, 0.2, 0.2, 0.2))
  
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
## Convert rownames from GeneSymbol to EnsID. 1 note.
addENSID <- function(dataframe, pc = TRUE) {
  if (pc) {
    annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only!  
  } else {
    annot <- geneInfo
  }
  
  # Match symbols in dataframe and annotation file
  sharedAnnotation <- annot[which(annot$Gene.Symbol%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$Gene.Symbol, rownames(dataframe))
  
  # Add ensID to a new column
  dataframe$ensID <- "-"
  dataframe$ensID[matches] <- sharedAnnotation$ensID # A sanity check was performed, which bound m$Approved.Symbol, and it always matched to rownames!
  dataframe <- dataframe[which(dataframe$ensID != "-"),]
  
  # Put ensID into rownames
  rownames(dataframe) <- dataframe$ensID
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "ensID")]
  return(dataframe)
}

load(file = "BrainCellularComposition/Preprocessed/geneInfo.rda")


countsData<-read.csv("sc_PUT_GLAND/GSE142653_pit_dev_5181_count.csv", header = TRUE, row.names = 1)

dat_ens = addENSID(countsData)

putgland <- CreateSeuratObject(counts = dat_ens, project = "GSE142653", min.cells = 3, min.features = 200)

metadata = read.csv('sc_PUT_GLAND/GSE142653_pit_dev_CellInfo.csv', row.names = 1, stringsAsFactors = T)

#select only cells that have metadata
putgland = putgland[,rownames(metadata)]

putgland = AddMetaData(putgland, metadata)

putgland$cellType = metadata$cell_type

putgland$cellID = colnames(putgland)

putgland = get_300_cells_per_celltype(putgland)
putgland = remove_rare_celltypes(putgland)

putgland <- NormalizeData(putgland, normalization.method = "LogNormalize", scale.factor = 10000)


sc = putgland
sc<- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType
DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0("put gland", ' (#: ', ncol(pbmc), ')' ))
#rename clusters
pbmc = RenameIdents(pbmc, "Thyrotrope" = "Thyrotrope/Lactotrope",  "Lactotrope" = "Thyrotrope/Lactotrope")

#sample some donors
pbmc  = equal_celltype_per_sample(pbmc)
pbmc$sampleID
saveRDS(pbmc, 'sc_PUT_GLAND/E-MTAB-5214_tabula_pituitarygland_seurat.rds')



pbmc$cellType = Idents(pbmc)
