library(plyr)
library(Seurat)


seurat = readRDS('/nfs/production/irene/ma/users/nnolte/DECONVOLUTION/tabula_sapiens/tabula_sapiens.rds')

print(seurat)


 
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
   "breast", "pancreas", "adipose tissue", "eye", "uterus", "minor salivary gland", 
   "skeletal muscle tissue", "atrium auricular region", "heart left ventricle", "eye", 
   "eye", "eye", "myometrium", "tibial artery", "tongue")


df_tissue = data.frame(seurat$tissue)
df_tissue$seurat.tissue <- mapvalues(df_tissue$seurat.tissue, from= tissue, to=ontology)

seurat$tissue = df_tissue$seurat.tissue

celltype = read.csv('files/cell_Type_list.tsv', sep = "\t", row.names = 1)

df_cellType= data.frame(seurat$cellType)
df_cellType$seurat.cellType <- mapvalues(df_cellType$seurat.cellType,  from= celltype$unique.seurat.cellType., to=celltype$modified_celltype)

seurat$cellType = df_cellType$seurat.cellType


saveRDS(seurat, '/nfs/production/irene/ma/users/nnolte/DECONVOLUTION/tabula_sapiens/tabula_sapiens_ontology_celltype_adapt.rds')


#write.csv(data.frame(unique(seurat$cellType)),'/nfs/production/irene/ma/users/nnolte/DECONVOLUTION/tabula_sapiens/cell_Type_list.tsv')