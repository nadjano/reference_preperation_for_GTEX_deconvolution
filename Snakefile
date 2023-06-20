
#Set config here
configfile: 'config.yaml'


#Helper functions to fetch inputs from other modules
def getMB(val):
        return int(val*1024)


rule all:
    input:
        # expand("Raw/Split/{source}.split.done", source=['TabulaSapiens', 'GTEx_sc','GutCellAtlas', 'IntegratedLungCellAtlas', 'HumanCellLandscapes',  'liver', 'ovary', 'Fallopian_tube'])
        # expand("UMAP/{tissue}_umap.png", tissue = ["HumanCellLandscapes_stomach", "HumanCellLandscapes_testis" ,"GTEx_sc_prostategland"])
        expand("Raw/Split/{source}.split.done", source = ["Ulrich_2021", "MacParland_2018", "Lengyel_2022", "GTEx_sc"]),
        # "UMAP/MacParland_2018_liver_umap.png"


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
        "Raw/Split/adult_human_20221007.loom.done"
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
        "Raw/Split/adult_human_20221007.loom.done"
    output:
        "Raw/Split/adult_human_20221007.split.done"
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
        "Raw/Split/{source}.split.done"
    conda:
        "env/loom.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])

    shell:
     """
     mkdir -p Raw/Split
     Rscript scripts/split_into_tissues.R {input} 
     touch {output}  
     """

rule reduce_cell_types:
    input:
        "Raw/Split/{tissue}_seurat.rds"
    output:
        "Raw/Split/{tissue}_seurat_curated.rds"
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
    log: "logs/{tissue}.log"
    input:
        "Raw/Split/{tissue}_seurat_curated.rds"
    output:
        "UMAP/{tissue}_umap.png"
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

 