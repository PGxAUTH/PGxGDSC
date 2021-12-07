# PGxGDSC
PGxGDSC is an automated pipeline for pharmacogenomic analysis of cancer drug response. Given a cancer type and therapeutic compound PGxGDSC identifies genes with mutually exclusive alterations in sensitive and resistance cancer cell models. Mutation profiles, gene expression and copy number variations are considered from COSMIC and CCLE repositories.

## Step-by-step description

### 1.Download prerequisite data to the working directory
Download the data files below and the R scripts found in **PGxGDSC/source/**. Then set the working directory to the downloaded data and R scripts.

**A. Drugs response data**
1. [GDSC1-dataset](http://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_25Feb20.xlsx)
2. [GDSC2-dataset](http://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx)
3. [Compounds-annotation](http://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/screened_compunds_rel_8.2.csv)
4. [List of all annotated models](https://cog.sanger.ac.uk/cmp/download/model_list_20210719.csv)

**B. COSMIC Cell Lines Project data files**
- You will need to login to download the files.
- COSMIC does not support older versions of the resource. As a result, *release v94* data files, which were used for our analysis, are not available for download. You can download the data files from COSMIC latest release, *release v95*, to use the methodology. However, the results of the analysis for the three use cases described in the article might be different from the ones provided in the **PGxGDSC/data/** file.

Go to the [COSMIC Cell Lines Project](https://cancer.sanger.ac.uk/cell_lines/download) and download the following files:
1. *Complete mutation data* > *Download Whole File* > *CosmicCLP_MutantExport.tsv.gz*
2. *Copy Number Data* > *Download Whole File* > *CosmicCLP_CompleteCNA.tsv.gz*
3. *Gene Expression* > *Download Whole File* > *CosmicCLP_CompleteGeneExpression.tsv.gz*

**C. DepMap/CCLE data files**

The CCLE data are available through DepMap and can be downloaded following the links below:
1. [CCLE_mutations](https://ndownloader.figshare.com/files/29125233)
2. [CCLE_expression_full](https://ndownloader.figshare.com/files/29124810)
3. [CCLE_gene_cn](https://ndownloader.figshare.com/files/29125230)

### 2.Run PGxGDSC/source/uploads_edits.R
The scripts loads the data files in the R workspace. This might take a while, so be patient.

### 3.Run PGxGDSC/source/main.R
Open the main.R script and configure the following parameters appropriately.
`drug_name <- `: Insert drug name as referenced in Genomics of Drug Sensitivity in Cancer database e.g. `drug_name <- "Afatinib"`.

**`dataset <- <dataset>`** 

*`<dataset>`: Insert Genomics of Drug Sensitivity in Cancer dataset of choice for the analysis. Available options are `dataset <- GDSC2` and `dataset <- GDSC1`.*

**`tissue <- <tissue>`**

*`<tissue>`: Insert TCGA classification for the tissue you are running the analysis e.g. `tissue <- "BRCA"` (BRCA is the TCGA classification for Breast invasive carcinoma)*

`getScript(x, y)`: Choose the type (*x*) and the source (*y*) of the data to be analyzed. Options for *x* are *"MUTATIONS"* (mutation data), *"EXPRESSION"* (gene expression data) and *"CN"* (copy number data). Options for *y* are *"COSMIC"* (COSMIC Cell Lines Project) and *"DEPMAP"* (DepMap/CCLE). e.g. `getScript("MUTATIONS", "COSMIC")`

If you want to run the analysis for another *drug_name*, *dataset* or *tissue*, you need first to run the command

`rm(list=setdiff(ls(), c("getScript", "GDSC1", "GDSC2", "compounds_annotation", "CCLE_exp", "CCLE_mutations", "CCLE_CN", "COSMIC_exp", "COSMIC_mutations", "COSMIC_CNVs")))`
then, insert a new *drug_name*, *dataset* or *tissue* and run `getScript(x, y)`.

Internet connection is required in order to download the data files and to run `source("KEGG.R")` and `source("pathway_interactions.R")`.

**Citation**: The work has been submitted for publication and is currently under review.




