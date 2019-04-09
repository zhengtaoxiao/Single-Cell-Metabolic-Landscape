Metabolic landscape of single cells in the tumor microenvironment
============
Introduction
------------
This pipeline is responsible for analyzing metabolic gene expression profile using single-cell RNA sequencing data. 

![pipeline](pipeline.png "Schematic representation of single-cell RNA-seq data analysis")

There are 7 main steps involved in data processing and analyses which organized in 7 different folders.

Requirements
------------
The required R packages can be installed using below command:

``` bash
Rscript install_requiredPackages.R 
```
Download and read the datasets
-----------------------------
``` bash
cd "1-ReadData"
bash download_dataset.sh
Rscript readData_head_neck.R
Rscript readData_melanoma.R
cd ../
```
The gene expression data and the annotation of the cell types would be stored as the R objects.

Imputation of the missig values
-------------------------------
``` bash
cd "2-Imputation"
Rscript impute_tpm.R melanoma 
Rscript geneNum_statistic.R melanoma
Rscript impute_tpm.R head_neck
Rscript geneNum_statistic.R head_neck
cd ../
```
This step uses the ["scImpute"](https://github.com/Vivianstats/scImpute) package to impute the missing values in gene expression profile. The gene number before and after imputation will be plotted as the boxplots.

Normalization and evaluation of different normalization methods 
---------------------------------------------------------------
``` bash
cd "3-Normalization"
Rscript normalization.R melanoma
Rscript normalization.R head_neck
cd ../
```
Four commonly used data normalization methods are applied on each dataset. The distribution of relative gene expression of each cell type will be ploted to evaluate and select the best normalization method.

Landscape of the metabolic gene expression profile
--------------------------------------------------
``` bash
cd "4-Clustering"
Rscript metabolic_landscape.R melanoma
Rscript metabolic_landscape.R head_neck
Rscript inter_tumor_distance.R melanoma
Rscript inter_tumor_distance.R head_neck
cd ../
```
The t-SNE algorithm will be performed in this step for visualizing metabolic gene expression in millions of cells. The spearman correlation matrix will aslo be generated to show the inter-tumor heterogeneity using metabolic genes.

Metabolic pathway activity
--------------------------
``` bash
cd 5-PathwayActivity
Rscript scRNA_pathway_activity.R melanoma
Rscript scRNA_pathway_activity.R head_neck
Rscript scRNA_pathway_activity_nonmetabolic.R melanoma
Rscript scRNA_pathway_activity_nonmetabolic.R head_neck
Rscript TCGA_pathway_activity.R
Rscript correlation_scRNA_TCGA.R
cd ..
```
This step will calculate the metabolic pathway activities for different cell types and for bulk tumor/normal samples using scRNA-seq dataset and TCGA dataset. The scatter plot for visualizing the correlation of pathway activities at single-cell resolution and bulk-tumor resolution will be generated. The violin plots will be generated to show the distribution of metabolic pathway activities in each cell types or bulk tumor/normal. To exclude that the imbalanced distribution of pathway activities in different cell types is due to artifactual, the distribution of non-metabolic pathway activities also be generated.

*The bulk RNA-seq data used here was downloaded from TCGA website, please see the instruction for data downloading and preprocessing in Data/TCGA/README.md* 

Metabolic pathway heterogeneity
-------------------------------
``` bash
cd 6-PathwayHeterogeneity
Rscript intra_malignant_variance_PCA.R melanoma
Rscript intra_malignant_variance_PCA.R head_neck
Rscript intra_non-malignant_variance_PCA.R melanoma
Rscript intra_non-malignant_variance_PCA.R head_neck
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R melanoma
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R head_neck
Rscript CCLE_OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R
Rscirpt Low_OXPHOS_glycolysis_hypoxia_signature.R melanoma
Rscript Low_OXPHOS_glycolysis_hypoxia_signature.R head_neck
cd ..
```
In this step, the PCA and GSEA analysis will be performed to investigate the pathway heterogeneity across cells in each tumor or cell type. The correlaton coefficients between the activities of OXPHOS, glycolysis and response to hypoxia will also be calculated using the scRNA gene expression data and the data from the cultured cell lines. The script "Low_OXPHOS_glycolysis_hypoxia_signature.R" is used for identifing the signature genes in OXPHOS/glycolysis/hypoxia low cells. The signature genes would be stored in the txt files,which can be used to perform the GO analysis on this website: http://metascape.org

Metabolic features of cell subtypes
-----------------------------------
``` bash
cd 7-MetabolicPhenotype
Rscript non-malignant_subtype.R melanoma
Rscript non-malignant_subtype.R head_neck
cd ..
```
The metabolic gene expression will be compared across T-cell subtypes and fibroblast subtypes using the GSEA analysis. 

Contact
-------
zhengtao.xiao@duke.edu