library(plyr)
library(Seurat)


seurat = readRDS('Raw/tabula_sapiens.rds')


source('scripts/R_functions.R')

seurat$cellType = seurat$cell_type

seurat$sampleID = seurat$donor_id

seurat$cellID = colnames(seurat)
 
tissue = c("coronary artery", "blood", "subcutaneous adipose tissue", "liver", "prostate gland", 
"spleen", "aorta", "uterus", "lung", "kidney", "large intestine", "lymph node", "endocrine pancreas",
 "thymus", "bladder organ", "inguinal lymph node", "small intestine", "trachea", "eye", "skin of abdomen",
  "skin of body", "skin of chest", "parotid gland", "bone marrow", "muscle tissue", "muscle of pelvic diaphragm", 
  "rectus abdominis muscle", "retinal neural layer", "anterior part of tongue", "posterior part of tongue", 
  "mammary gland", "exocrine pancreas", "adipose tissue", "lacrimal gland", "endometrium", "submandibular gland",
   "muscle of abdomen", "cardiac atrium", 
"cardiac ventricle", "conjunctiva", "sclera", "cornea", "myometrium", "vasculature", "tongue")


ontology = c("coronary artery", "blood", "subcutaneous adipose tissue", "liver", "prostate gland", "spleen", 
"aorta", "uterus", "lung", "cortex of kidney", "sigmoid colon", "lymph node", "pancreas", "thymus", 
"urinary bladder", "lymph node", "small intestine Peyer's patch", "trachea", "eye", "lower leg skin",
 "lower leg skin", "lower leg skin", "minor salivary gland", "bone marrow", "skeletal muscle tissue",
  "muscle of pelvic diaphragm", "rectus abdominis muscle", "eye", "tongue", "tongue",
   "breast", "pancreas", "subcutaneous adipose tissue", "eye", "uterus", "minor salivary gland", 
   "skeletal muscle tissue", "atrium auricular region", "heart left ventricle", "eye", 
   "eye", "eye", "uterus", "tibial artery", "tongue")


df_tissue = data.frame(seurat$tissue)
df_tissue$seurat.tissue <- mapvalues(df_tissue$seurat.tissue, from= tissue, to=ontology)

seurat$tissue = df_tissue$seurat.tissue

celltype = read.csv('tabels/cell_Type_list.tsv', sep = "\t", row.names = 1)

df_cellType= data.frame(seurat$cellType)
df_cellType$seurat.cellType <- mapvalues(df_cellType$seurat.cellType,  from= celltype$unique.seurat.cellType., to=celltype$modified_celltype)

seurat$cellType = df_cellType$seurat.cellType



GTEx_tissues = c("coronary artery", "blood", "subcutaneous adipose tissue", "prostate gland", "spleen", 
"aorta", "uterus", "cortex of kidney", "pancreas", "urinary bladder",  "minor salivary gland",
 "pancreas", "subcutaneous adipose tissue", "uterus", "minor salivary gland", 
 "atrium auricular region", "tibial artery", "lower leg skin")

for (tissue in unique(GTEx_tissues)){


      output = gsub(" ", "", paste0('Processed/', tissue, '_seurat.rds'))
      split_sc_per_tissue(seurat, tissue, output)
      
}




#saveRDS(seurat, 'Raw/tabula_sapiens_ontology_celltype_adapt.rds')


#write.csv(data.frame(unique(seurat$cellType)),'/nfs/production/irene/ma/users/nnolte/DECONVOLUTION/tabula_sapiens/cell_Type_list.tsv')