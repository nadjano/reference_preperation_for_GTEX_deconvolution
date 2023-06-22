



remove_rare_celltypes = function(seurat){
    # this function removes cells from cell types with less than 100 cells
    # from SeuratObject
    row_count <- as.data.frame(x = table(seurat$cellType))
    non_rare = row_count[row_count$Freq >100 ,]
    seurat <- subset(seurat, subset = cellType %in% paste(non_rare$Var1))
    return(seurat)
  
}



get_300_cells_per_celltype <- function(seurat){
     ## downsamples cell to max 300 cells per celltype
  
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
  
  # Add ensID to a new column
  dataframe$ensID <- "-"
  dataframe$ensID[matches] <- sharedAnnotation$ensID # A sanity check was performed, which bound m$Approved.Symbol, and it always matched to rownames!
  dataframe <- dataframe[which(dataframe$ensID != "-"),]
  
  # Put ensID into rownames
  rownames(dataframe) <- dataframe$ensID
  dataframe <- dataframe[,-which(colnames(dataframe) %in% "ensID")]
  return(dataframe)
}


equal_celltype_per_sample <- function(seurat){
  
  donor = c('TSP1', 'TSP2', 'TSP3', 'TSP4', 'TSP5')
  seurat$sampleID = sample(donor, size = length(colnames(seurat)), replace=T, prob =c(0.2, 0.2, 0.2, 0.2, 0.2))
  
  return(seurat)
}


split_sc_per_tissue <- function(seurat, tissue, output){
    
    
    if (!file.exists(output)) {
        #seurat = readRDS(seurat)

        cellSplits <- SplitObject(seurat, split.by = "tissue")
                
            
        for (i in 1:length(levels(seurat$tissue))){
        
        
            sc = cellSplits[[i]]
            sc$tissue = droplevels(x = sc$tissue)
            sc$cellType = droplevels(x = sc$cellType)
            this_tissue = levels(unique((sc$tissue)))
            if (this_tissue ==  tissue){
                sc = get_300_cells_per_celltype(sc)
                #sc = equal_celltype_per_sample(sc)
                sc = remove_rare_celltypes(sc)
                saveRDS(sc, output)

            }
         
    }

    } else {
                print(paste(output, "already exists. Skipping..."))
            }
}


get_semantic_tag <- function(tissue, ontology) {
  tissue <- gsub('_| ', '+', tissue)
  
  url <- "www.ebi.ac.uk/spot/zooma/v2/api/services/annotate?propertyValue="
  res <- GET(paste0(url, tissue) )
  stop_for_status(res)
  data <- fromJSON(rawToChar(res$content))
  #filter to only get organism part ontologies
  if (length(data) != 0){
    data = data[ grepl(ontology,(data$semanticTags) ),]
    #check confidence of zooma maping
    if (nrow(data) > 0){
      data = data[1, ]
      if (data$confidence != 'LOW'){
        semanticTag <- basename(data$semanticTags[[1]])
        return(semanticTag)
      } else {
        stop(paste('No high/good confidence mapping found for', tissue))
      }} else {
        paste('No mapping found for', tissue)
      }
    
  } else {
    paste('No mapping found for', tissue)
  }
  
}

