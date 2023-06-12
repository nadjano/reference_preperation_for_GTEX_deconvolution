library(stringr)
library(dplyr)
library(plyr)
library(ggplot2)
library(Seurat)


load('Raw/SRA826293_SRS4181130.sparse.RData')

cluster = read.csv('Raw/SRA826293_SRS4181130.clusters.txt', header = T, sep = ",")
celltypes = read.csv('Raw/inferredcellType.tsv', sep = "\t")

str_split_fixed(str_split_fixed(row, "ENSG",)[,-1], ".", -1)

merged = inner_join(cluster, celltypes,by  = 'clusterID' , )


first = str_split_fixed(rownames(seurart), "ENSG", 2)[,2]

ens_rownames = paste0('ENSG', str_split_fixed(first,fixed("."), 2)[,1])


rownames(sm) = ens_rownames

seurat = CreateSeuratObject(sm)

seurat = seurat[, merged$Barcode]

seurat$cellType = merged$Inferred.cell.type

#remove unknown celltypes
seurat = subset(seurat, subset = cellType != 'Unknown')

#remove celltypes with less than 100 cells
seurat = remove_rare_celltypes(seurat)

#downsample to 300 cells per celltype
seurat= get_300_cells_per_celltype(seurat)

sc = seurat
sc<- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType
DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0("testis", ' (#: ', ncol(seurat), ')' ))

pbmc$cellID = colnames(pbmc)

#sample some sampleIDs as we dont have this information
pbmc = equal_celltype_per_sample(pbmc)

saveRDS(pbmc, 'Processed/testis_seurat.rds')
