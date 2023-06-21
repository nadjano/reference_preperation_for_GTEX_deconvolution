
#Set config here
configfile: 'config.yaml'


#Helper functions to fetch inputs from other modules
def getMB(val):
        return int(val*1024)


wildcard_constraints:
    source="^(?!.*Linnarsson_2022).*"
    
rule all:
    input:
        # expand("Raw/Split/{source}.split.done", source=['TabulaSapiens', 'GTEx_sc','GutCellAtlas', 'IntegratedLungCellAtlas', #'HumanCellLandscapes',  'liver', 'ovary', 'Fallopian_tube'])
        # expand("UMAP/{tissue}_umap.png", tissue = ["HumanCellLandscapes_stomach", "HumanCellLandscapes_testis" ,"GTEx_sc_prostategland"])
        #expand("Raw/Split/{source}.split.done", source = ["Ulrich_2021", "MacParland_2018", "Lengyel_2022", "GTEx_sc"]),
        'UMAP/MacParland_2018_liver_umap.png',
        'UMAP/TabulaSapiens_coronary-artery_umap.png',
        'UMAP/TabulaSapiens_blood_umap.png',
        'UMAP/TabulaSapiens_subcutaneous-adipose-tissue_umap.png',
        'UMAP/TabulaSapiens_prostate-gland_umap.png',
        'UMAP/TabulaSapiens_aorta_umap.png',
        'UMAP/TabulaSapiens_uterus_umap.png',
        'UMAP/TabulaSapiens_kidney_umap.png',
        'UMAP/TabulaSapiens_endocrine-pancreas_umap.png',
        'UMAP/TabulaSapiens_bladder-organ_umap.png',
        'UMAP/TabulaSapiens_parotid-gland_umap.png',
        'UMAP/TabulaSapiens_mammary-gland_umap.png',
        'UMAP/TabulaSapiens_muscle-tissue_umap.png',
        'UMAP/TabulaSapiens_exocrine-pancreas_umap.png',
        'UMAP/TabulaSapiens_adipose-tissue_umap.png',
        'UMAP/TabulaSapiens_endometrium_umap.png',
        'UMAP/TabulaSapiens_cardiac-atrium_umap.png',
        'UMAP/TabulaSapiens_cardiac-ventricle_umap.png',
        'UMAP/TabulaSapiens_myometrium_umap.png',
        'UMAP/GTEx_sc_esophagus-muscularis-mucosa_umap.png',
        'UMAP/GTEx_sc_mucosa_umap.png',
        'UMAP/GutCellAtlas_large-intestine_umap.png',
        'UMAP/GutCellAtlas_small-intestine_umap.png',
        'UMAP/Wiedemann_2023_skin-of-body_umap.png',
        'UMAP/Ulrich_2021_fallopian-tube_umap.png',
        'UMAP/Lengyel_2022_ovary_umap.png',
        'UMAP/HumanCellLandscapes_stomach_umap.png',
        'UMAP/IntegratedLungCellAtlas_lungparenchyma_umap.png',
        'UMAP/HumanCellLandscapes_testis_umap.png',
        'UMAP/HumanCellLandscapes_thyroid-gland_umap.png'
        #'UMAP/Linnarsson_2022_Pons_umap.png'

rule adrenal_gland_reference:
    input:
        "Raw/hAG_cnts_ensIDs_seurat.rds"

    output:
        "Raw/Split/adrenal_gland_seurat.rds"
    
    conda:
        "env.yaml"

    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    
    shell:
        """
        Rscript scripts/sc_FALLOPIANTUBE.R
        """

rule split_brain:
    input:
        loom="Raw/adult_human_20221007.loom",
        clusters="Raw/cluster_annotation.xlsx"
    output:
        "Raw/Split/adult_human_20221007.split.done"
    conda:
        "env/loom_py.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     mkdir -p Raw/Split
     python scripts/devide_loom.py {input.loom} {input.clusters}  
     touch Raw/Split/adult_human_20221007.loom.done
     """

rule prepare_brain:
    input:
        "Raw/Split/adult_human_20221007.split.done"
    output:
        expand("Raw/Split/Linnarsson_2022_{brain_tissue}_seurat.rds", brain_tissue = ["Cerebellum", "Cerebralcortex", "Cerebralnuclei", "Hippocampus", "Midbrain", "Pons", "Spinalcord", "Thalamus", "Hypothalamus"])
    conda:
        "env/loom.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     mkdir -p Raw/Split
     Rscript scripts/loom.R  
     """

rule split_into_tissues:
    input:
        "Raw/{source}.rds"
    output:
        "Raw/Split/{source}_{tissue}_seurat.rds"
    conda:
        "env/loom.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])

    shell:
     """
     mkdir -p Raw/Split
     Rscript scripts/split_into_tissues.R {input} {wildcards.tissue}
     """

rule reduce_cell_types:
    input:
        "Raw/Split/{source}_{tissue}_seurat.rds"
    output:
        "Raw/Split/{source}_{tissue}_seurat_curated.rds"
    conda:
        "env/scONTO.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     Rscript scripts/reduceCellTypes.R {input}   
     """

rule reduce_cell_types_brain:
    input:
        "Raw/Split/Linnarsson_2022_{brain_tissue}_seurat.rds"
    output:
        "Raw/Split/Linnarsson_2022_{brain_tissue}_seurat_curated.rds"
    conda:
        "env/scONTO.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     Rscript scripts/reduceCellTypes.R {input}   
     """


rule UMAP_plots:
    """
    Rule to plot UMAP plots for quality control of references. Created two UMAP plots
    in one png file with old and reduced cell type labels.
    """
    conda: "env/UMAP.yaml"
    log: "logs/{source}_{tissue}.log"
    input:
        "Raw/Split/{source}_{tissue}_seurat_curated.rds"
    output:
        "UMAP/{source}_{tissue}_umap.png"
    resources:  mem_mb=getMB(config['mem_gb']['DWLS'])
    shell: 
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        mkdir -p UMAP
        Rscript {workflow.basedir}/scripts/makeUMAPplots.R {input} {output}
        """

rule UMAP_plots_brain:
    """
    Rule to plot UMAP plots for quality control of references. Created two UMAP plots
    in one png file with old and reduced cell type labels.
    """
    conda: "env/UMAP.yaml"
    log: "logs/{brain_tissue}.log"
    input:
        "Raw/Split/Linnarsson_2022_{brain_tissue}_seurat_curated.rds"
    output:
        "UMAP/Linnarsson_2022_{brain_tissue}_umap.png"
    resources:  mem_mb=getMB(config['mem_gb']['DWLS'])
    shell: 
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        mkdir -p UMAP
        Rscript {workflow.basedir}/scripts/makeUMAPplots.R {input} {output}
        """

# rule rename_final_files:

# shell:
# """
# mkdir -p FinalOutput
# mv Fallopian_tube_UBERON_0016632_seurat_curated.rds FinalOutput/GTEx_v8-Fallopian_Tube_seurat.rds
# mv GTEx_sc_prostategland_seurat_curated.rds FinalOutput/GTEx_v8-Prostate_seurat.rds
# mv TabulaSapiens_muscletissue_seurat_curated.rds GTEx_v5-Muscle_-_Skeletal_seurat.rds
# mv TabulaSapiens_adiosetissue_seurat_curated.rds  GTEx_v5-Adipose_-_Subcutaneous_seurat.rds
# GTEx_v5-Adipose_-_Subcutaneous_seurat.rds GTEx_v5-Adipose_-_Visceral_Omentum
# """

 