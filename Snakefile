
#Set config here
configfile: 'config.yaml'


#Helper functions to fetch inputs from other modules
def getMB(val):
        return int(val*1024)

ruleorder: split_brain > prepare_brain > split_into_tissues > reduce_cell_types_brain > reduce_cell_types > UMAP_plots_brain > UMAP_plots

# rule all:
#     input:
#         'UMAP/MacParland_2018_liver_umap.png',
#         'UMAP/TabulaSapiens_coronary-artery_umap.png',
#         'UMAP/TabulaSapiens_blood_umap.png',
#         'UMAP/TabulaSapiens_spleen_umap.png',
#         'UMAP/TabulaSapiens_subcutaneous-adipose-tissue_umap.png',
#         'UMAP/TabulaSapiens_aorta_umap.png',
#         'UMAP/TabulaSapiens_uterus_umap.png',
#         'UMAP/TabulaSapiens_kidney_umap.png',
#         'UMAP/TabulaSapiens_bladder-organ_umap.png',
#         'UMAP/TabulaSapiens_parotid-gland_umap.png',
#         'UMAP/TabulaSapiens_mammary-gland_umap.png',
#         'UMAP/TabulaSapiens_muscle-tissue_umap.png',
#         'UMAP/TabulaSapiens_adipose-tissue_umap.png',
#         'UMAP/TabulaSapiens_endometrium_umap.png',
#         'UMAP/TabulaSapiens_cardiac-atrium_umap.png',
#         'UMAP/TabulaSapiens_cardiac-ventricle_umap.png',
#         'UMAP/TabulaSapiens_myometrium_umap.png',
#         'UMAP/GTEx_sc_esophagus-muscularis-mucosa_umap.png',
#         'UMAP/GTEx_sc_mucosa_umap.png',
#         'UMAP/GutCellAtlas_large-intestine_umap.png',
#         'UMAP/GutCellAtlas_small-intestine_umap.png',
#         'UMAP/Wiedemann_2023_skin-of-body_umap.png',
#         'UMAP/Ulrich_2021_fallopian-tube_umap.png',
#         'UMAP/Lengyel_2022_ovary_umap.png',
#         'UMAP/HumanCellLandscapes_stomach_umap.png',
#         'UMAP/HumanCellLandscapes_pancreas_umap.png',
#         'UMAP/HumanCellLandscapes_adrenal-gland_umap.png',
#         'UMAP/IntegratedLungCellAtlas_lung-parenchyma_umap.png',
#         'UMAP/HumanCellLandscapes_testis_umap.png',
#         'UMAP/Linnarsson_2022_Cerebellum_umap.png',
#         'UMAP/Linnarsson_2022_Cerebralnuclei_umap.png',
         # 'UMAP/Linnarsson_2022_Cerebralcortex_umap.png',
#         'UMAP/Linnarsson_2022_Midbrain_umap.png',
#         'UMAP/Linnarsson_2022_Hippocampus_umap.png',
#         'UMAP/Linnarsson_2022_Spinalcord_umap.png',
#         'UMAP/Linnarsson_2022_Thalamus_umap.png',
#         'UMAP/Linnarsson_2022_Hypothalamus_umap.png'

rule all:
    input:
        'all.done'


rule adrenal_gland_reference:
    input:
        "Raw/hAG_cnts_ensIDs_seurat.rds"

    output:
        "Raw/adrenal_gland_seurat.rds"
    
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
        "Raw/adult_human_20221007.split.done"
    conda:
        "env/loom_py.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     mkdir -p Split
     python scripts/devide_loom.py {input.loom} {input.clusters}  
     touch Split/adult_human_20221007.loom.done
     """

rule prepare_brain:
    input:
        "Split/adult_human_20221007.split.done"
    output:
        expand("Split/Linnarsson_2022_{brain_tissue}_seurat.rds", brain_tissue = ["Cerebellum", "Cerebralcortex", "Cerebralnuclei", "Hippocampus", "Midbrain", "Pons", "Spinalcord", "Thalamus", "Hypothalamus"])
    conda:
        "env/loom.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     mkdir -p Split
     Rscript scripts/loom.R  
     """

rule split_into_tissues:
    input:
        "Raw/{source}.rds"
    output:
        "Split/{source}_{tissue}_seurat.rds"
    conda:
        "env/loom.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])

    shell:
     """
     mkdir -p Raw
     Rscript scripts/split_into_tissues.R {input} {wildcards.tissue}
     """
     
rule reduce_cell_types_brain:
    input:
        "Split/Linnarsson_2022_{brain_tissue}_seurat.rds"
    output:
        "Split/Linnarsson_2022_{brain_tissue}_seurat_curated.rds"
    conda:
        "env/scONTO.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     Rscript scripts/reduceCellTypes.R {input}   
     """

rule reduce_cell_types:
    input:
        "Split/{source}_{tissue}_seurat.rds"
    output:
        "Split/{source}_{tissue}_seurat_curated.rds"
    conda:
        "env/scONTO.yaml"
    resources:
        mem_mb=getMB(config['mem_gb']['DWLS'])
    shell:
     """
     Rscript scripts/reduceCellTypes.R {input}   
     """

rule UMAP_plots_brain:
    """
    Rule to plot UMAP plots for quality control of references. Created two UMAP plots
    in one png file with old and reduced cell type labels.
    """
    conda: "env/UMAP.yaml"
    log: "logs/{brain_tissue}.log"
    input:
        "Split/Linnarsson_2022_{brain_tissue}_seurat_curated.rds"
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

rule UMAP_plots:
    """
    Rule to plot UMAP plots for quality control of references. Created two UMAP plots
    in one png file with old and reduced cell type labels.
    """
    conda: "env/UMAP.yaml"
    log: "logs/{source}_{tissue}.log"
    input:
        "Split/{source}_{tissue}_seurat_curated.rds"
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

rule rename_final_files:
    input:
        'UMAP/MacParland_2018_liver_umap.png',
        'UMAP/TabulaSapiens_coronary-artery_umap.png',
        'UMAP/TabulaSapiens_blood_umap.png',
        'UMAP/TabulaSapiens_spleen_umap.png',
        'UMAP/TabulaSapiens_subcutaneous-adipose-tissue_umap.png',
        'UMAP/TabulaSapiens_aorta_umap.png',
        'UMAP/TabulaSapiens_uterus_umap.png',
        'UMAP/TabulaSapiens_kidney_umap.png',
        'UMAP/TabulaSapiens_bladder-organ_umap.png',
        'UMAP/TabulaSapiens_parotid-gland_umap.png',
        'UMAP/TabulaSapiens_mammary-gland_umap.png',
        'UMAP/TabulaSapiens_muscle-tissue_umap.png',
        'UMAP/TabulaSapiens_adipose-tissue_umap.png',
        'UMAP/TabulaSapiens_endometrium_umap.png',
        'UMAP/TabulaSapiens_cardiac-atrium_umap.png',
        'UMAP/TabulaSapiens_cardiac-ventricle_umap.png',
        'UMAP/TabulaSapiens_myometrium_umap.png',
        'UMAP/GTEx_sc_esophagus-muscularis-mucosa_umap.png',
        'UMAP/GTEx_sc_mucosa_umap.png',
        'UMAP/GTEx_sc_prostate-gland_umap.png',
        'UMAP/GutCellAtlas_large-intestine_umap.png',
        'UMAP/GutCellAtlas_small-intestine_umap.png',
        'UMAP/Wiedemann_2023_skin-of-body_umap.png',
        'UMAP/Ulrich_2021_fallopian-tube_umap.png',
        'UMAP/Lengyel_2022_ovary_umap.png',
        'UMAP/HumanCellLandscapes_stomach_umap.png',
        'UMAP/HumanCellLandscapes_pancreas_umap.png',
        'UMAP/HumanCellLandscapes_adrenal-gland_umap.png',
        'UMAP/IntegratedLungCellAtlas_lung-parenchyma_umap.png',
        'UMAP/HumanCellLandscapes_testis_umap.png',
        'UMAP/Linnarsson_2022_Cerebellum_umap.png',
        'UMAP/Linnarsson_2022_Cerebralnuclei_umap.png',
        'UMAP/Linnarsson_2022_Midbrain_umap.png',
        'UMAP/Linnarsson_2022_Hippocampus_umap.png',
        'UMAP/Linnarsson_2022_Spinalcord_umap.png',
        'UMAP/Linnarsson_2022_Thalamus_umap.png',
        'UMAP/Linnarsson_2022_Hypothalamus_umap.png',
        'UMAP/Linnarsson_2022_Cerebralcortex_umap.png'
    output: 'all.done'
    params: 
        sampleName='GTEx_v8'
    shell:
        """
        mkdir -p FinalOutput

        cp Split/GTEx_sc_esophagus-muscularis-mucosa_seurat_curated.rds FinalOutput/{params.sampleName}-Esophagus_-_Muscularis_seurat.rds
        cp Split/GTEx_sc_esophagus-muscularis-mucosa_seurat_curated.rds FinalOutput/{params.sampleName}-Esophagus_-_Gastroesophageal_Junction_seurat.rds
        cp Split/GTEx_sc_mucosa_seurat_curated.rds FinalOutput/{params.sampleName}-Esophagus_-_Mucosa_seurat.rds
        cp Split/GTEx_sc_prostate-gland_seurat_curated.rds FinalOutput/{params.sampleName}-Prostate_seurat.rds
        cp Split/HumanCellLandscapes_thyroid-gland_seurat_curated.rds FinalOutput/{params.sampleName}-Thyroid_seurat.rds
        cp Split/IntegratedLungCellAtlas_lung-parenchyma_seurat_curated.rds FinalOutput/{params.sampleName}-Lung_seurat.rds
        cp Split/Lengyel_2022_ovary_seurat_curated.rds FinalOutput/{params.sampleName}-Ovary_seurat.rds
        cp Split/MacParland_2018_liver_seurat_curated.rds FinalOutput/{params.sampleName}-Liver_seurat.rds
        cp Split/TabulaSapiens_adipose-tissue_seurat_curated.rds FinalOutput/{params.sampleName}-Adipose_-_Visceral_Omentum_seurat.rds
        cp Split/TabulaSapiens_aorta_seurat_curated.rds FinalOutput/{params.sampleName}-Artery_-_Aorta_seurat.rds
        cp Split/TabulaSapiens_bladder-organ_seurat_curated.rds FinalOutput/{params.sampleName}-Bladder_seurat.rds
        cp Split/TabulaSapiens_blood_seurat_curated.rds FinalOutput/{params.sampleName}-Whole_Blood_seurat.rds
        cp Split/TabulaSapiens_blood_seurat_curated.rds FinalOutput/{params.sampleName}-Cells_-_EBV-transformed_lymphocytes_seurat.rds
        cp Split/TabulaSapiens_cardiac-atrium_seurat_curated.rds FinalOutput/{params.sampleName}-Heart_-_Atrial_Appendage_seurat.rds
        cp Split/TabulaSapiens_cardiac-ventricle_seurat_curated.rds FinalOutput/{params.sampleName}-Heart_-_Left_Ventricle_seurat.rds
        cp Split/TabulaSapiens_coronary-artery_seurat_curated.rds FinalOutput/{params.sampleName}-Artery_-_Coronary_seurat.rds
        cp Split/TabulaSapiens_coronary-artery_seurat_curated.rds FinalOutput/{params.sampleName}-Artery_-_Tibial_seurat.rds
        cp Split/TabulaSapiens_endometrium_seurat_curated.rds FinalOutput/{params.sampleName}-Cervix_-_Endocervix_seurat.rds
        cp Split/TabulaSapiens_endometrium_seurat_curated.rds FinalOutput/{params.sampleName}-Vagina_seurat.rds
        cp Split/TabulaSapiens_endometrium_seurat_curated.rds FinalOutput/{params.sampleName}-Cervix_-_Ectocervix_seurat.rds
        cp Split/TabulaSapiens_kidney_seurat_curated.rds FinalOutput/{params.sampleName}-Kidney_-_Cortex_seurat.rds
        cp Split/TabulaSapiens_kidney_seurat_curated.rds FinalOutput/{params.sampleName}-Kidney_-_Medulla_seurat.rds
        cp Split/TabulaSapiens_mammary-gland_seurat_curated.rds FinalOutput/{params.sampleName}-Breast_-_Mammary_Tissue_seurat.rds
        cp Split/TabulaSapiens_muscle-tissue_seurat_curated.rds FinalOutput/{params.sampleName}-Muscle_-_Skeletal_seurat.rds
        cp Split/TabulaSapiens_parotid-gland_seurat_curated.rds FinalOutput/{params.sampleName}-Minor_Salivary_Gland_seurat.rds
        cp Split/TabulaSapiens_subcutaneous-adipose-tissue_seurat_curated.rds FinalOutput/{params.sampleName}-Adipose_-_Subcutaneous_seurat.rds
        cp Split/TabulaSapiens_uterus_seurat_curated.rds FinalOutput/{params.sampleName}-Uterus_seurat.rds
        cp Split/Ulrich_2021_fallopian-tube_seurat_curated.rds FinalOutput/{params.sampleName}-Fallopian_Tube_seurat.rds
        cp Split/Wiedemann_2023_skin-of-body_seurat_curated.rds FinalOutput/{params.sampleName}-Skin_-_Not_Sun_Exposed_Suprapubic_seurat.rds
        cp Split/Wiedemann_2023_skin-of-body_seurat_curated.rds FinalOutput/{params.sampleName}-Skin_-_Sun_Exposed_Lower_leg_seurat.rds
        cp Split/Wiedemann_2023_skin-of-body_seurat_curated.rds FinalOutput/{params.sampleName}-Cells_-_Cultured_fibroblasts_seurat.rds
        cp Split/TabulaSapiens_spleen_seurat_curated.rds FinalOutput/{params.sampleName}-Spleen_seurat.rds
        cp Split/HumanCellLandscapes_stomach_seurat_curated.rds FinalOutput/{params.sampleName}-Stomach_seurat.rds
        cp Split/HumanCellLandscapes_pancreas_seurat_curated.rds FinalOutput/{params.sampleName}-Pancreas_seurat.rds
        cp Split/HumanCellLandscapes_testis_seurat_curated.rds FinalOutput/{params.sampleName}-Testis_seurat.rds
        cp Split/HumanCellLandscapes_adrenal-gland_seurat_curated.rds FinalOutput/{params.sampleName}-Adrenal_Gland_seurat.rds
        cp Split/GutCellAtlas_small-intestine_seurat_curated.rds FinalOutput/{params.sampleName}-Small_Intestine_-_Terminal_Ileum_seurat.rds
        cp Split/GutCellAtlas_large-intestine_seurat_curated.rds FinalOutput/{params.sampleName}-Colon_-_Sigmoid_seurat.rds
        cp Split/GutCellAtlas_large-intestine_seurat_curated.rds FinalOutput/{params.sampleName}-Colon_-_Transverse_seurat.rds
        cp Split/Linnarsson_2022_Hippocampus_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Amygdala_seurat.rds
        cp Split/Linnarsson_2022_Cerebralcortex_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Anterior_cingulate_cortex_BA24_seurat.rds
        cp Split/Linnarsson_2022_Cerebralnuclei_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Caudate_basal_ganglia_seurat.rds
        cp Split/Linnarsson_2022_Cerebellum_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Cerebellar_Hemisphere_seurat.rds
        cp Split/Linnarsson_2022_Cerebellum_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Cerebellum_seurat.rds
        cp Split/Linnarsson_2022_Cerebralcortex_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Cortex_seurat.rds
        cp Split/Linnarsson_2022_Cerebralcortex_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Frontal_Cortex_BA9_seurat.rds
        cp Split/Linnarsson_2022_Hippocampus_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Hippocampus_seurat.rds
        cp Split/Linnarsson_2022_Hypothalamus_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Hypothalamus_seurat.rds
        cp Split/Linnarsson_2022_Cerebralnuclei_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Nucleus_accumbens_basal_ganglia_seurat.rds
        cp Split/Linnarsson_2022_Cerebralnuclei_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Putamen_basal_ganglia_seurat.rds
        cp Split/Linnarsson_2022_Spinalcord_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Spinal_cord_cervical_c-1_seurat.rds
        cp Split/Linnarsson_2022_Midbrain_seurat_curated.rds FinalOutput/{params.sampleName}-Brain_-_Substantia_nigra_seurat.rds
        cp Split/Linnarsson_2022_Spinalcord_seurat_curated.rds FinalOutput/{params.sampleName}-Nerve_-_Tibial_seurat.rds
        touch all.done
        """