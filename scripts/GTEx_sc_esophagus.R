library(plyr)
#library(scOntoMatch)
library(Seurat)
library(ggplot2)
library(stringr)

####
# select esophagusmuscularis and esophagusmucosa from GTEx_sc 
# rename celltypes
# remove rare cell types
# sample 300 cells per celltype

source('scripts/R_functions.R')

GTEx_sc = readRDS('Raw/GTEx_sc.rds')
print(head(colnames(GTEx_sc)))
print(head(rownames(GTEx_sc)))
print(head(GTEx_sc@meta.data))


GTEx_sc$donorID = GTEx_sc$donor_id
GTEx_sc$cellID = colnames(GTEx_sc)
GTEx_sc$cellType = GTEx_sc$cell_type


GTEx_sc$tissue = factor(str_split(colnames(GTEx_sc), '-', 2, simplify = T)[, 2])
tissues = unique(GTEx_sc$tissue)
print(tissues)

for (tissue in tissues){
  
  split_sc_per_tissue(GTEx_sc, tissue, gsub(" ", "", paste0('Processed/', tissue, '_seurat.rds')))

}





# GTEx_sc = subset(GTEx_sc, subset = tissue %in% c('esophagusmucosa', 'esophagusmuscularis'))
# df = data.frame(GTEx_sc$`Broad cell type`)
# unique(GTEx_sc$`Granular cell type`)
# #GTEx_sc = GTEx_sc[,GTEx_sc$tissue == 'heart']

# org_cellType = c("Adipocyte", "Endothelial cell (lymphatic)", "Endothelial cell (vascular)", "Epithelial cell (Hillock)", "Epithelial cell (alveolar type I)", "Epithelial cell (alveolar type II)", "Epithelial cell (basal keratinocyte)", "Epithelial cell (basal)", "Epithelial cell (ciliated)", "Epithelial cell (club)", "Epithelial cell (cornified keratinocyte)", "Epithelial cell (luminal)", "Epithelial cell (mature keratinocyte)", "Epithelial cell (squamous)", "Epithelial cell (suprabasal keratinocyte)", "Epithelial cell (suprabasal)", "Fibroblast", "ICCs", "Immune (B cell)", "Immune (DC)", "Immune (DC/macrophage)", "Immune (Langerhans)", "Immune (NK cell)", "Immune (T cell)", "Immune (alveolar macrophage)", "Immune (mast cell)", "Melanocyte", "Mucous cell", "Myocyte (NMJ-rich)", "Myocyte (cardiac)", "Myocyte (cardiac, cytoplasmic)", "Myocyte (sk. muscle)", "Myocyte (sk. muscle, cytoplasmic)", "Myocyte (smooth muscle)", "Myoepithelial (basal)", "Myofibroblast", "Neuroendocrine", "Neuronal", "Pericyte/SMC", "Satellite cell", "Schwann cell", "Sebaceous gland cell", "Sweat gland cell", "Unknown")
# new_cellType = c("fat cell", "endothelial cell", "endothelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "epithelial cell", "fibroblast", "interstitial cell of Cajal", "B cell", "dendritic cell", "macrophage", "Langerhans cell", "mature NK T cell", "T cell", "macrophage", "mast cell", "melanocyte", "epithelial cell", "muscle cell", "cardiac muscle cell", "cardiac muscle cell", "muscle cell", "muscle cell", "smooth muscle cell", "basal cell", "smooth muscle cell", "neuroendocrine cell", "neuron", "pericyte cell", "skeletal muscle satellite stem cell", "Schwann cell", "sebaceous gland cell", "sweat secreting cell", "NA")

# # get ontology cellType names
# GTEx_sc$ontology_cellType =  mapvalues((df$GTEx_sc..Broad.cell.type.), from = (org_cellType), to = new_cellType)

# GTEx_sc$cellType = GTEx_sc$ontology_cellType

# #sample per celltype 300 cells
# GTEx_sc$tissue = droplevels(x = GTEx_sc$tissue)
# GTEx_sc$cellType = droplevels(x = GTEx_sc$cellType)

# GTEx_sc = get_300_cells_per_celltype(GTEx_sc)
# GTEx_sc = remove_rare_celltypes(GTEx_sc)

# sc = GTEx_sc
# sc<- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(sc)
# pbmc <- ScaleData(sc, features = all.genes)

# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# pbmc <- RunUMAP(pbmc, dims = 1:10)

# Idents(pbmc) = pbmc$cellType
# DimPlot(pbmc, label = TRUE)   + NoLegend() +ggtitle(paste0('esophagusmuscularis', ' (#: ', ncol(pbmc), ')' ))

# pbmc$sampleID = pbmc$individual
# pbmc$cellID = colnames(pbmc)

# saveRDS(pbmc, 'Processed/esophagusmucosa_seurat.rds')

# #eso = readRDS( 'GTEx_sc/E-MTAB-5214_tabula_esophagusmucosa_seurat.rds')
# #rownames(eso)



