library(tidyverse)

getScript <- function(x, y) {
  
  if(x == "MUTATIONS" & y == "COSMIC") {source("MUTATIONS_COSMIC.R")}
  if(x == "MUTATIONS" & y == "DEPMAP") {source("MUTATIONS_DEPMAP.R")}
  if(x == "EXPRESSION" & y == "COSMIC") {source("EXPRESSION_COSMIC.R")}
  if(x == "EXPRESSION" & y == "DEPMAP") {source("EXPRESSION_DEPMAP.R")}
  if(x == "CN" & y == "COSMIC") {source("CN_COSMIC.R")}
  if(x == "CN" & y == "DEPMAP") {source("CN_DEPMAP.R")}
  
}

drug_name <- "Afatinib"
dataset <- GDSC2
tissue <- "BRCA"
core_geneset <- list()

getScript("MUTATIONS", "COSMIC")

#### Clear all to Search for a new drug_name / dataset / tissue
rm(list=setdiff(ls(), c("getScript", "GDSC1", "GDSC2", "compounds_annotation", "CCLE_exp", "CCLE_mutations", "CCLE_CN", "COSMIC_exp", "COSMIC_mutations", "COSMIC_CNVs")))

#### Sum up results
core_geneset <- bind_rows(core_geneset)
core_geneset

#### Search for drug's target pathways
source("KEGG.R")
kegg_TARGET_PATHWAY

#### Search for pathway interactions
results <- list()
source("pathway_interactions.R")
results <- bind_rows(results)
