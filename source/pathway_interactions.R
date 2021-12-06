library(tidyverse)

library(clusterProfiler)

kegg_TARGET_PATHWAY_GENES_SYMBOL <- kegg_TARGET_PATHWAY_GENES

for(i in 1:length(kegg_TARGET_PATHWAY_GENES_SYMBOL)){
  kegg_TARGET_PATHWAY_GENES_SYMBOL[[i]] <- bitr(kegg_TARGET_PATHWAY_GENES_SYMBOL[[i]], "ENTREZID", "SYMBOL", "org.Hs.eg.db")[,2]
}

pathway_genes <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- core_geneset$`Genes` %in% x[[i]]
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

pathway_gene <- pathway_genes(kegg_TARGET_PATHWAY_GENES_SYMBOL)
names(pathway_gene) <- kegg_TARGET_PATHWAY

for(i in 1:length(pathway_gene)){
  pathway_gene[[i]] <- replace(pathway_gene[[i]], pathway_gene[[i]]==TRUE, names(pathway_gene)[[i]])
}


core_geneset_filter <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- core_geneset %>% mutate("Pathway" = x[[i]])
    listofx[[i]] <- df
    
  }
  return(listofx)
}

core_geneset_filter <- core_geneset_filter(pathway_gene)
core_geneset_filter <- bind_rows(core_geneset_filter)
core_geneset_filter <- core_geneset_filter %>% filter(Pathway!=FALSE)
core_geneset_filter <- core_geneset_filter %>% mutate("Pathway target"="Pathway gene", "PPI score"=NA, "Core/Extended network"="Core")

results$"Pathway genes" <- core_geneset_filter

###############################################

library("STRINGdb")
string_db_1 <- STRINGdb$new( version="11.5", species=9606, score_threshold=900, input_directory="")
string_db_2 <- STRINGdb$new( version="11.5", species=9606, score_threshold=950, input_directory="")


core_geneset_string <- string_db_1$map( core_geneset, "Genes", removeUnmappedRows = TRUE )


TARGET_PATHWAY_GENES <- kegg_TARGET_PATHWAY_GENES

for(i in 1:length(TARGET_PATHWAY_GENES)){
  TARGET_PATHWAY_GENES[[i]] <- as.data.frame(TARGET_PATHWAY_GENES[[i]])
  colnames(TARGET_PATHWAY_GENES[[i]]) <- "genes"
}

for(i in 1:length(TARGET_PATHWAY_GENES)){
  TARGET_PATHWAY_GENES[[i]] <- string_db_1$map( TARGET_PATHWAY_GENES[[i]], "genes", removeUnmappedRows = TRUE )
}

################## DIRECT ################## 

neighbors <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- string_db_1$get_neighbors(x[[i]])
    listofx[[i]] <- vector
    
  }
  return(listofx)
}


number_of_neighbors <- vector()
for(i in 1:length(core_geneset_string$STRING_id)){
  number_of_neighbors[[i]] <- length(string_db_1$get_neighbors(core_geneset_string$STRING_id[[i]]))
}

core_geneset_string <- core_geneset_string %>% mutate(number_of_neighbors=number_of_neighbors)
core_geneset_string <- core_geneset_string %>% filter(number_of_neighbors!=0)
neighbors_list <- neighbors(core_geneset_string$STRING_id)
names(neighbors_list) <- core_geneset_string$STRING_id


neighborsfrompathways <- function(x, y) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- filter(y, y$STRING_id %in% x[[i]] == TRUE)
    listofx[[i]] <- df
    
  }
  return(listofx)
}

listofinteractions <- list()

for(i in 1:length(TARGET_PATHWAY_GENES)){
  listofinteractions[[i]] <- neighborsfrompathways(neighbors_list, TARGET_PATHWAY_GENES[[i]] )
}


names(listofinteractions) <- kegg_TARGET_PATHWAY
StringIDs <- core_geneset_string$STRING_id


for(i in 1:length(listofinteractions)){
  names(listofinteractions[[i]]) <- StringIDs
}


for(i in 1:length(listofinteractions)){
  listofinteractions[[i]] <- listofinteractions[[i]] %>% keep(~ nrow(.x) != 0)
}

listofinteractions <- compact(listofinteractions)

if (length(listofinteractions)!=0) {

interactions_score_1 <- function(x, y) {
  listofx <- vector()
  for(i in 1:length(x)){
    vector <- string_db_1$get_interactions( c(y, x[[i]]  ))[1,3]
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

interactions_score_2 <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- x[[i]] %>% mutate(combined_score = interactions_score_1(x[[i]]$STRING_id, names(x)[[i]]))
    listofx[[i]] <- df
    
  }
  return(listofx)
}

interactions_score_3 <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    list <- interactions_score_2(x[[i]])
    listofx[[i]] <- list
    
  }
  return(listofx)
}

direct_interactions_list <- interactions_score_3(listofinteractions)

for(i in 1:length(direct_interactions_list)){
  names(direct_interactions_list[[i]]) <- names(listofinteractions[[i]])
}

names(direct_interactions_list) <- names(listofinteractions)


mutate_protein_name <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- x[[i]] %>% mutate("Protein" = names(x)[[i]])
    listofx[[i]] <- df
    
  }
  return(listofx)
}

for(i in 1:length(direct_interactions_list)){
  direct_interactions_list[[i]] <- mutate_protein_name(direct_interactions_list[[i]])
}

for(i in 1:length(direct_interactions_list)){
  direct_interactions_list[[i]] <- bind_rows(direct_interactions_list[[i]])
}

mutate_pathway_name <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- x[[i]] %>% mutate("Pathway" = names(x)[[i]])
    listofx[[i]] <- df
    
  }
  return(listofx)
}

direct_interactions_list <- mutate_pathway_name(direct_interactions_list)
direct_interactions_list <- bind_rows(direct_interactions_list)

genes_entrezid_symbols <- bitr(direct_interactions_list$genes, "ENTREZID", "SYMBOL", "org.Hs.eg.db")
colnames(genes_entrezid_symbols) <- c("genes", "Pathway target")

direct_interactions_list <- left_join(direct_interactions_list, genes_entrezid_symbols)
direct_interactions_list[,1] <- direct_interactions_list$Pathway
direct_interactions_list[,2] <- direct_interactions_list$`Pathway target`
direct_interactions_list <- direct_interactions_list %>% transmute("Pathway"=genes, "Pathway target"=STRING_id, "PPI score"=combined_score, "STRING_id"=Protein)

"Direct interactions" <- right_join(core_geneset_string, direct_interactions_list, by="STRING_id")
`Direct interactions` <- `Direct interactions` %>% transmute("Genes"=Genes, "Cell Lines Class"=Cell.Lines.Class, "Source"=Source, "Type of analysis"=Type.of.analysis,  "Type of alteration"=Type.of.alteration, "Pathway"=Pathway, "Pathway target"=`Pathway target`, "PPI score"=`PPI score`, "Core/Extended network"="Core")

`Direct interactions` <- unique(`Direct interactions`)
results$"Direct interactions" <- `Direct interactions`
}

################## INDIRECT ################## 

for(i in 1:length(neighbors_list)){
  neighbors_list[[i]] <- as.data.frame(neighbors_list[[i]])
  colnames(neighbors_list[[i]]) <- "STRING_id"
}


mutate_function <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- x[[i]] %>% mutate(core = names(x)[[i]])
    listofx[[i]] <- df
    
  }
  return(listofx)
}


neighbors_list <- mutate_function(neighbors_list)
names(neighbors_list) <- StringIDs

extended_geneset <- bind_rows(neighbors_list)

core_geneset_string_selected <- core_geneset_string[,c("Genes", "STRING_id")]
colnames(core_geneset_string_selected) <- c("gene", "core")
extended_geneset <- left_join(extended_geneset, core_geneset_string_selected)


neighbors2 <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- string_db_2$get_neighbors(x[[i]])
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

number_of_neighbors2 <- vector()
for(i in 1:length(extended_geneset$STRING_id)){
  number_of_neighbors2[[i]] <- length(string_db_2$get_neighbors(extended_geneset$STRING_id[[i]]))
}

extended_geneset <- extended_geneset %>% mutate(number_of_neighbors=number_of_neighbors2)
extended_geneset <- extended_geneset %>% filter(number_of_neighbors!=0)
neighbors_list2 <- neighbors2(extended_geneset$STRING_id)
names(neighbors_list2) <- extended_geneset$STRING_id


neighborsfrompathways2 <- function(x, y) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- filter(y, y$STRING_id %in% x[[i]] == TRUE)
    listofx[[i]] <- df
    
  }
  return(listofx)
}

listofinteractions2 <- list()

for(i in 1:length(TARGET_PATHWAY_GENES)){
  listofinteractions2[[i]] <- neighborsfrompathways2(neighbors_list2, TARGET_PATHWAY_GENES[[i]] )
}

names(listofinteractions2) <- kegg_TARGET_PATHWAY
StringIDs2 <- extended_geneset$STRING_id

for(i in 1:length(listofinteractions2)){
  names(listofinteractions2[[i]]) <- StringIDs2
}

for(i in 1:length(listofinteractions2)){
  listofinteractions2[[i]] <- listofinteractions2[[i]] %>% keep(~ nrow(.x) != 0)
}

listofinteractions2 <- compact(listofinteractions2)

if (length(listofinteractions2)!=0) {

interactions_score_1_2 <- function(x, y) {
  listofx <- vector()
  for(i in 1:length(x)){
    vector <- string_db_2$get_interactions( c(y, x[[i]]  ))[1,3]
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

interactions_score_2_2 <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    df <- x[[i]] %>% mutate(combined_score = interactions_score_1_2(x[[i]]$STRING_id, names(x)[[i]]))
    listofx[[i]] <- df
    
  }
  return(listofx)
}

interactions_score_3_2 <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    list <- interactions_score_2_2(x[[i]])
    listofx[[i]] <- list
    
  }
  return(listofx)
}

indirect_interactions_list <- interactions_score_3_2(listofinteractions2)

for(i in 1:length(indirect_interactions_list)){
  names(indirect_interactions_list[[i]]) <- names(listofinteractions2[[i]])
}

names(indirect_interactions_list) <- names(listofinteractions2)


for(i in 1:length(indirect_interactions_list)){
  indirect_interactions_list[[i]] <- mutate_protein_name(indirect_interactions_list[[i]])
}

for(i in 1:length(indirect_interactions_list)){
  indirect_interactions_list[[i]] <- bind_rows(indirect_interactions_list[[i]])
}

indirect_interactions_list <- mutate_pathway_name(indirect_interactions_list)
indirect_interactions_list <- bind_rows(indirect_interactions_list)

genes_entrezid_symbols <- bitr(indirect_interactions_list$genes, "ENTREZID", "SYMBOL", "org.Hs.eg.db")
colnames(genes_entrezid_symbols) <- c("genes", "Pathway target")

indirect_interactions_list <- left_join(indirect_interactions_list, genes_entrezid_symbols)
indirect_interactions_list[,1] <- indirect_interactions_list$Pathway
indirect_interactions_list[,2] <- indirect_interactions_list$`Pathway target`
indirect_interactions_list <- indirect_interactions_list %>% transmute("Pathway"=genes, "Pathway target"=STRING_id, "PPI score"=combined_score, "STRING_id"=Protein)

indirect_interactions_list <- unique(indirect_interactions_list)

indirect_interactions_list <- left_join(indirect_interactions_list, extended_geneset)

library("biomaRt")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
indirect_interactions_list$STRING_id <- as.data.frame(str_split_fixed(indirect_interactions_list$STRING_id, "9606.", 2))$V2
extended_geneset_names <-  getBM(attributes = c('hgnc_symbol', 'ensembl_peptide_id'), filters = 'ensembl_peptide_id',  values = indirect_interactions_list$STRING_id, mart = ensembl)
colnames(extended_geneset_names)[2] <- "STRING_id"

missing_gene_names <- as.data.frame(as.data.frame(str_split_fixed(extended_geneset$STRING_id, "9606.", 2))$V2)
colnames(missing_gene_names) <- "STRING_id"
missing_gene_names <- anti_join(missing_gene_names, extended_geneset_names)
missing_gene_names <- missing_gene_names %>% transmute('hgnc_symbol'=STRING_id, "STRING_id"=STRING_id)
extended_geneset_names <- rbind(extended_geneset_names, missing_gene_names)

indirect_interactions_list <- left_join(indirect_interactions_list, extended_geneset_names)
indirect_interactions_list <- indirect_interactions_list %>% transmute("Pathway"=Pathway, `Pathway target`=`Pathway target`, `PPI score`=`PPI score`, "STRING_id"=core, "Core/Extended network"=hgnc_symbol)
core_geneset_string <- core_geneset_string %>% transmute("Genes"=Genes, "Cell Lines Class"=Cell.Lines.Class, "Source"=Source, "Type of analysis"=Type.of.analysis,  "Type of alteration"=Type.of.alteration, "STRING_id"=STRING_id)

"Indirect interactions" <- right_join(core_geneset_string, indirect_interactions_list, by="STRING_id")
`Indirect interactions` <- `Indirect interactions`[,-6]

`Indirect interactions` <- unique(`Indirect interactions`)
results$"Indirect interactions" <- `Indirect interactions`
}


