
#Set config here
configfile: 'config.yaml'


#Helper functions to fetch inputs from other modules
def getMB(val):
        return int(val*1024)


# rule ovary_reference:
#     input:
#         "Raw/ovary.rds"

#     output:
#         "Processed/ovary_seurat.rds"
    
#     conda:
#         "env.yaml"
    
#     resources:
#         mem_mb=getMB(config['mem_gb']['DWLS'])
    
#     shell:
#         """
#         Rscript scripts/sc_OVARY.R
#         """


# rule fallopiantube_reference:
#     input:
#         "Raw/FallopianTube.rds"

#     output:
#         "Processed/fallopiantube_seurat.rds"
    
#     conda:
#         "env.yaml"

#     resources:
#         mem_mb=getMB(config['mem_gb']['DWLS'])
    
#     shell:
#         """
#         Rscript scripts/sc_FALLOPIANTUBE.R
#         """


rule split_into_tissues:
    input:
        "Raw/GTEx_sc_esophagus.rds"
    output:
        "Raw/Split/GTEx_sc.split.done"
    conda:
        "env.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])

    script:
     """
     mkdir -p Raw/Split
     Rscript scripts/split_into_tissues.R {input} 
     touch {output}  
     """

rule reduce_cell_types:
    input:
        "Raw/Split/liver_seurat.rds"
    output:
        "Raw/Split/liver_seurat_curated.rds"
    conda:
        "env/scONTO.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     Rscript scripts/reduceCellTypes.R {input}   
     """