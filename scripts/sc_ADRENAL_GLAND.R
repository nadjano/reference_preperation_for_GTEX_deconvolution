#filename <- 'sc_ADRENAL_GLAND/hAG_cnts_ensIDs.h5ad'



#convert h5ad to seurat
#suppressMessages(library(devtools))
#devtools::install_github('cellgeni/sceasy')
#suppressMessages(library(sceasy))
#suppressMessages(library(Seurat))
#sceasy::convertFormat(filename, from = "anndata", to = "seurat", outFile = paste0(sub(".h5ad", "", filename), "_seurat.rds"))


library('Seurat')
library(dplyr)
library(plyr)
library(ggplot2)

source('scripts/R_functions.R')

seurat = readRDS('Raw/hAG_cnts_ensIDs_seurat.rds')

#gene_seurat = readRDS('sc_ADRENAL_GLAND/hAG_cnts_seurat.rds')
seurat$samples
old_label = c('Progenitor_1', 'Macrophages_2', 'Zfasciculata_3','Chromaffin_4', 'Zglomerulosa_5', 'Endothellial_6','Mesenchymal_7','Zreticularis_8','Cortex_9','Tcells_10')    
new_label = c('progenitor', 'macrophages', 'zona fasciculata','chromaffin', 'zona glomerulosa ', 'endothellial','mesenchymal','zona reticularis','cortex','T.cells')  

seurat$cellType <- mapvalues(seurat$PAGODA_hc,   from= old_label, to=new_label)

seurat$sampleID = seurat$samples
seurat$cellID = colnames(seurat)

seurat = remove_rare_celltypes(seurat)
seurat = get_300_cells_per_celltype(seurat)

seurat<- NormalizeData(seurat)

sc<- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(sc)
pbmc <- ScaleData(sc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


pbmc <- RunUMAP(pbmc, dims = 1:10)

Idents(pbmc) = pbmc$cellType
#DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0("adrenal gland", ' (#: ', ncol(pbmc), ')' ))
#rename clusters
saveRDS(pbmc, 'Processed/adrenalgland_seurat.rds')
