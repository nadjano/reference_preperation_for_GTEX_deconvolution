
#Set config here
configfile: 'config.yaml'


#Helper functions to fetch inputs from other modules
def getMB(val):
        return int(val*1024)





rule ovary_reference:
    input:
        "Raw/ovary.rds"

    output:
        "Processed/ovary_seurat.rds"
    
    conda:
        "env.yaml"
    
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    
    shell:
        """
        Rscript scripts/sc_OVARY.R
        """


rule fallopiantubereference_reference:
    input:
        "Raw/FallopianTube.rds"

    output:
        "Processed/fallopiantube_seurat.rds"
    
    conda:
        "env.yaml"

    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    
    shell:
        """
        Rscript scripts/sc_FALLOPIANTUBE.R
        """