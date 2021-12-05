library(tidyverse)
library(httr)
library(readxl)

# Upload GDSC datasets
GET("ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_25Feb20.xlsx", write_disk(tf <- tempfile(fileext = ".xlsx")))
GDSC1 <- read_xlsx(tf)
GET("ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_25Feb20.xlsx", write_disk(tf <- tempfile(fileext = ".xlsx")))
GDSC2 <- read_xlsx(tf)
compounds_annotation <- read.csv("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/screened_compunds_rel_8.2.csv")

# Upload sample info from Sanger's Cell Model Passports
sample_info <- read.csv("model_list_20210719.csv") # https://cog.sanger.ac.uk/cmp/download/model_list_20210719.csv
sample_info <- sample_info %>% select(BROAD_ID, model_id)

#### MUTATIONS DATA ####
# Upload complete mutation data from COSMIC Cell Lines Project
COSMIC_mutations <- read_tsv("CosmicCLP_MutantExport.tsv") # https://cancer.sanger.ac.uk/cell_lines/download
## Remove ENSTs
COSMIC_mutations <- COSMIC_mutations %>% mutate(Gene_name_edited = as.data.frame(str_split_fixed(COSMIC_mutations$`Gene name`, "_ENST", 2))$V1)
## Mutations with FATHMM score > 0.5 are considered PATHOGENIC
COSMIC_mutations$`FATHMM prediction` <- ifelse(COSMIC_mutations$`FATHMM score`>0.5 & COSMIC_mutations$`FATHMM score`<0.7, "PATHOGENIC", COSMIC_mutations$`FATHMM prediction`)

# Upload mutation data from DepMap
CCLE_mutations <- read.csv("https://ndownloader.figshare.com/files/29125233")
## Add sample info
colnames(sample_info) <- c("DepMap_ID", "Sanger_Model_ID")
CCLE_mutations <- left_join(CCLE_mutations, sample_info, by="DepMap_ID")

#### EXPRESSION DATA ####
# Upload Gene Expression data from COSMIC Cell Lines Project
COSMIC_exp <- read_tsv("CosmicCLP_CompleteGeneExpression.tsv") # https://cancer.sanger.ac.uk/cell_lines/download
# Edit 
Cell_Lines <- as.factor(COSMIC_exp$SAMPLE_NAME)
Cell_Lines <- levels(Cell_Lines)
getlist <- function(Cell_Lines) {
  listofcell_lines <- list()
  for(i in 1:length(Cell_Lines)){
    df <- COSMIC_exp %>% filter(SAMPLE_NAME == Cell_Lines[i])
    listofcell_lines[[i]] <- df
    
  }
  return(listofcell_lines)
}
CellLinesList <- getlist(Cell_Lines = Cell_Lines)
CellLinesList <- lapply(CellLinesList, function(x) x[(names(x) %in% c("GENE_NAME", "Z_SCORE"))])
COSMIC_exp <- CellLinesList %>% reduce(left_join, by = "GENE_NAME")
colnames(COSMIC_exp)[2:971] <- Cell_Lines

# Upload expression data from DepMap
CCLE_exp <- read.csv("CCLE_expression_full.csv") # "https://ndownloader.figshare.com/files/29124810"
# Edit
colnames(sample_info) <- c("X", "model_id")
CCLE_exp <- left_join(CCLE_exp, sample_info)
CCLE_exp[,1] <- CCLE_exp$model_id
CCLE_exp <- CCLE_exp %>% select(-model_id)
CCLE_exp <- as.data.frame(t(CCLE_exp))
Cell_Lines <- as.character(CCLE_exp[1,])
colnames(CCLE_exp) <- Cell_Lines
CCLE_exp <- CCLE_exp[-1,]
CCLE_exp[,1:1377] <- apply(CCLE_exp[,1:1377], 2, function(x) as.numeric(as.character(x)))

# Upload Copy Number Data from COSMIC Cell Lines Project
COSMIC_CNVs <- read_tsv("CosmicCLP_CompleteCNA.tsv") # https://cancer.sanger.ac.uk/cell_lines/download
# Edit
COSMIC_CNVs <- COSMIC_CNVs %>% mutate(Gene_name_edited = as.data.frame(str_split_fixed(COSMIC_CNVs$gene_name, "_ENST", 2))$V1)
COSMIC_CNVs <- COSMIC_CNVs %>% mutate(Gene_MUT_TYPE = paste(COSMIC_CNVs$Gene_name_edited, COSMIC_CNVs$MUT_TYPE, sep=", "))
COSMIC_CNVs <- COSMIC_CNVs %>% filter(!str_detect(`Chromosome:G_Start..G_Stop`, "Y:"))
COSMIC_CNVs <- COSMIC_CNVs %>% filter(!str_detect(`Chromosome:G_Start..G_Stop`, "X:"))

# Upload Copy Number Data from DepMap
CCLE_CN <- read.csv("CCLE_gene_cn.csv") # https://ndownloader.figshare.com/files/29125230
# Edit
colnames(sample_info) <- c("X", "model_id")
CCLE_CN <- left_join(CCLE_CN, sample_info)
CCLE_CN[,1] <- CCLE_CN$model_id
CCLE_CN <- CCLE_CN %>% select(-model_id)
CCLE_CN <- as.data.frame(t(CCLE_CN))
Cell_Lines <- as.character(CCLE_CN[1,])
colnames(CCLE_CN) <- Cell_Lines
CCLE_CN <- CCLE_CN[-1,]
CCLE_CN[,1:1741] <- apply(CCLE_CN[,1:1741], 2, function(x) as.numeric(as.character(x)))


