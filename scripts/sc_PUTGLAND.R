
library(Seurat)
library(ggplot2)
library("plyr")

source('scripts/R_functions.R')

#seurat = Read10X('sc_PUT_GLAND/')


## Convert rownames from GeneSymbol to EnsID. 1 note.


load(file = "tabels/geneInfo.rda")


countsData<-read.csv("Raw/GSE142653_pit_dev_5181_count.csv", header = TRUE, row.names = 1)

dat_ens = addENSID(countsData)

putgland <- CreateSeuratObject(counts = dat_ens, project = "GSE142653", min.cells = 3, min.features = 200)

metadata = read.csv('Raw/GSE142653_pit_dev_CellInfo.csv', row.names = 1, stringsAsFactors = T)

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

saveRDS(pbmc, 'Processed/putgland_seurat.rds')



