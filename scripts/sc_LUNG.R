library(Seurat)
library(ggplot2)
## Convert rownames from GeneSymbol to EnsID. 1 note.
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
addENSID <- function(dataframe, pc = TRUE) {
  if (pc) {
    annot <- geneInfo[which(geneInfo$Biotype == "protein_coding"),] # Note 1: selects protein-coding genes only!  
  } else {
    annot <- geneInfo
  }
  
  # Match symbols in dataframe and annotation file
  sharedAnnotation <- annot[which(annot$Gene.Symbol%in%rownames(dataframe)),]
  matches <- match(sharedAnnotation$Gene.Symbol, rownames(dataframe))
  print(matches)
  # Add ensID to a new column
  dataframe$ensID <- "-"
  dataframe$ensID[matches] <- sharedAnnotation$ensID # A sanity check was performed, which bound m$Approved.Symbol, and it always matched to rownames!
  print(dataframe)
  #dataframe <- dataframe[which(dataframe$ensID != "-"),]
  
  # Put ensID into rownames
  #rownames(dataframe) <- dataframe$ensID
 
  #dataframe <- dataframe[,-which(colnames(dataframe) %in% "ensID")]
  return(dataframe)
  
}
genes = read.csv('sc_LUNG/features.tsv', row.names = 1, header = FALSE)

load(file = "BrainCellularComposition/Preprocessed/geneInfo.rda")
genes_ens = addENSID(genes)
df = genes_ens
df$ensID[df$ensID == '-'] <- rownames(df)

#write.table(genes_ens, sep = '\t', row.names = T, col.names = T, file = 'sc_LUNG/ens_features.tsv')

counts = ReadMtx('sc_LUNG/gene_sorted-lung_expression_data.mtx', cells = 'sc_LUNG/barcodes.tsv', features = 'sc_LUNG/features.tsv',  feature.column = 1)


rownames(counts) = df$ensID



seurat = CreateSeuratObject(counts)

metadata = read.csv('sc_LUNG/metadata.tsv', sep = '\t', row.names = 1)

seurat = AddMetaData(seurat, metadata)



seurat_healthy =subset(seurat, subset = disease__ontology_label == 'normal')
seurat_healthy$cellType = seurat_healthy$cell_type_main

seurat_healthy_300  = get_300_cells_per_celltype(seurat_healthy)

sc<- FindVariableFeatures(seurat_healthy_300, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType


DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle('lung')

saveRDS(pbmc, 'sc_LUNG/lung_300.rds')
#plot(dimplot)