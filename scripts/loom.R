#library(scater)
library(SeuratDisk)
#library(reticulate)
library(Seurat)
#library(loomR)
library(dplyr)
library(plyr)
library(ggplot2)
library(Matrix)
library(stringr)

source('scripts/R_functions.R')

# load loom file that was split into brain regions and downsamples to 
# 5000 cells per region


files = list.files('Raw/brain_split/', pattern = '*.loom')
regions = str_remove(files, '_5000.loom')
print(regions)
#tissues = c("lung", "blood", "pancreas", "liver", 'coronaryartery', 'cortexofkidney', 'EBV-transformedlymphocyte', 'endocervix', 'ectocervix', 'uterus')


for (region in regions){
  
  seurat = LoadLoom(paste0('Raw/brain_split/', region, '_5000.loom' ))


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

#DimPlot(pbmc, label = TRUE)  +ggtitle('Cerebellum_8000')

# rename myofibroblast to fibroblast
#pbmc$cellType = mapvalues(pbmc$cellType, from = 'myofibroblast cell', to = 'fibroblast')
  pbmc$cellID= pbmc$sampleID
  pbmc$sampleID = colnames(pbmc)

# Save the Seurat object to a file
  saveRDS(pbmc, file = paste0("Processed/", region, '_seurat.rds'))


}