# Reference Preparation for GTEx Deconvolution

This repository contains a Snakemake pipeline to create a single-cell reference for deconvolving all GTEx tissues.

## Prerequisites
Before running the pipeline, make sure you have the following:

- Snakemake installed on your system.
- Clone this repository to your local machine.
- Create a directory named `Raw` and download the following single-cell references into it:

### Tabula Sapiens

- Download all cell data from [here](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) and save it as `Raw/TabulaSapiens.rds`.

### Atlas of the Developing Human Brain

- Download the brain cell atlas from [this link](https://storage.googleapis.com/linnarsson-lab-human/adult_human_20221007.loom).
- Download the cluster annotations from [this link](https://github.com/linnarsson-lab/adult-human-brain/raw/main/tables/cluster_annotation.xlsx).

### Integrated Human Lung Cell Atlas

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293) and save it as `Raw/IntegratedLungCellAtlas.rds`.

### Ovary

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/d36ca85c-3e8b-444c-ba3e-a645040c6185) and save it as `Raw/Lengyel_2022.rds`.

### Fallopian Tube

- Download the scRNA-seq analysis of healthy human fallopian tubes .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/fc77d2ae-247d-44d7-aa24-3f4859254c2c) and save it as `Raw/Ulrich_2021.rds`.

### Liver

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/bd5230f4-cd76-4d35-9ee5-89b3e7475659) and save it as `Raw/MacParland_2018.rds`.

### Gut Cell Atlas

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/e33ffcd3-7cbf-4b8c-b0f4-85587ad5019a) and save it as `Raw/GutCellAtlas.rds`.
- Download the Pediatric .rds fiel CellxGene [here](https://cellxgene.cziscience.com/collections/17481d16-ee44-49e5-bcf0-28c0780d8c4a) and save is at `Raw/Elmentaite_2020.rds`

### Pituitary Gland

- Download the data from [GSE142653](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142653) and [geneInfo.rda](https://github.com/Voineagulab/BrainCellularComposition/blob/main/Files/geneInfo.rda?raw=true).
- Obtain the files `GSE142653_pit_dev_5181_count.csv` and `GSE142653_pit_dev_CellInfo.csv`.

### Skin

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/e/da684768-fb01-455b-9f0f-b63a3e2f844f.cxg/) and save it as `Raw/Wiedemann_2023.rds`.

### Construction of a Human Cell

 Landscape at Single-Cell Level

- Download the .rds file from CellxGene [here](https://cellxgene.cziscience.com/collections/38833785-fac5-48fd-944a-0f62a4c23ed1) and save it as `Raw/HumanCellLandscapes.rds`.

## Additional Files
- Download the cell type ontology file from [here](http://purl.obolibrary.org/obo/cl/cl-basic.obo).

## Running the Pipeline
To run the pipeline, use the following command:
```
snakemake --profile lsf
```

Please make sure you have met all the prerequisites and have the required files in the appropriate directories before running the pipeline.
