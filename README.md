# reference_preperation_for_GTEX_deconvolution

This rep contains the script to make single cell reference to deconvolve GTEx tissue by tissue. 
The scripts are in scripts.


Create a directory 'Raw' and download the following references:

## single cell references

### Tabula Sapiens

+ downlaoaded all cell from here (https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5)
and save as `Raw/TabulaSapiens.rds`

### Atlas of the developing human brain:

* [Brain cell Atlas](https://storage.googleapis.com/linnarsson-lab-human/adult_human_20221007.loom)
* [cluster annotations](https://github.com/linnarsson-lab/adult-human-brain/raw/main/tables/cluster_annotation.xlsx)


### Integrated Human Lung Cell Atlas

+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)
and save as `Raw/IntegratedLungCellAtlas.rds`

### Ovary

+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/collections/d36ca85c-3e8b-444c-ba3e-a645040c6185)
and save as `Raw/Lengyel_2022.rds`

### Fallopian Tube

+ [Downlad scRNA-seq analysis of healthy human fallopian tubes .rds file from CellxGene](https://cellxgene.cziscience.com/collections/fc77d2ae-247d-44d7-aa24-3f4859254c2c) and save as `Raw/Ulrich_2021.rds`

### Liver

+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/collections/bd5230f4-cd76-4d35-9ee5-89b3e7475659) and save as `MacParland_2018.rds`


### Gut Cell Atlas
+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/collections/e33ffcd3-7cbf-4b8c-b0f4-85587ad5019a)
and save as `GutCellAtlas.rds`

### Putiary Gland

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142653 and https://github.com/Voineagulab/BrainCellularComposition/blob/main/Files/geneInfo.rda?raw=true (to tabels )to change gene names to ENS IDs

get GSE142653_pit_dev_5181_count.csv and GSE142653_pit_dev_CellInfo.csv

### Skin
+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/e/da684768-fb01-455b-9f0f-b63a3e2f844f.cxg/)
and save as `Wiedemann_2023.rds`

### Construction of a human cell landscape at single-cell level
+ [Downlad .rds file from CellxGene](https://cellxgene.cziscience.com/collections/38833785-fac5-48fd-944a-0f62a4c23ed1)
and save as `HumanCellLandscapes.rds`
