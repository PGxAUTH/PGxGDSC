library(tidyverse)
library(httr)

library("KEGGREST")
library("EnrichmentBrowser")


compounds_distinct <- compounds_annotation %>% distinct(DRUG_NAME, SYNONYMS)
compounds_distinct <- compounds_distinct %>% filter(is.na(SYNONYMS)==FALSE)
compounds_distinct <- filter(compounds_distinct, SYNONYMS!="")
compounds_distinct <- distinct(compounds_distinct, DRUG_NAME, .keep_all=TRUE)
# Berzosertib (VE-822)
compounds_distinct[272, 2] <- "Berzosertib"


compounds_distinct_drug_names <- compounds_distinct$DRUG_NAME
compounds_distinct_synonyms <- compounds_distinct$SYNONYMS
getSynonyms <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- unlist(strsplit(x[[i]], split=", "))
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

Synonyms <- getSynonyms(x=compounds_distinct_synonyms)
names(Synonyms) <- compounds_distinct_drug_names


getKeggFindResults <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- keggFind("drug", x[[i]])
    listofx[[i]] <- vector
    
  }
  return(listofx)
}


keggFindResult <- keggFind("drug", drug_name)

{
  if ((length(keggFindResult)==0) & (drug_name %in% compounds_distinct_drug_names==FALSE)) {stop("Drug not found.")}
  
}

ifelse(length(keggFindResult)==0, ifelse(length(compact(getKeggFindResults(x=Synonyms[[drug_name]])))!=0, keggFindResult <- compact(getKeggFindResults(x=Synonyms[[drug_name]]))[[1]], keggFindResult <- vector()), keggFindResult)

{
  if (length(keggFindResult)==0) {stop("Drug not found.")}
  
}

keggFindResult <- rownames(as.data.frame(keggFindResult))[1]

checkAntineoplastic <- lapply(keggGet(keggFindResult)[[1]], function(ch) grep("Antineoplastic", ch))
checkAntineoplastic <- compact(checkAntineoplastic)

{
  if (length(checkAntineoplastic)==0) {stop("Drug not found.")}
  
}

library(berryFunctions)

{
  if (is.error(keggGet(keggFindResult)[[1]]$TARGET$PATHWAY)==TRUE) {stop("Target pathways not found.")}
  
}

kegg_TARGET_PATHWAY <- keggGet(keggFindResult)[[1]]$TARGET$PATHWAY
kegg_TARGET_PATHWAY <- as.data.frame(str_split_fixed(kegg_TARGET_PATHWAY, "\\(", 2))$V1

gethsagenes <- function(x) {
  listofx <- list()
  for(i in 1:length(x)){
    vector <- keggGet(x[[i]])[[1]]$GENE[dput(seq(1,length(keggGet(x[[i]])[[1]]$GENE),2))]
    listofx[[i]] <- vector
    
  }
  return(listofx)
}

kegg_TARGET_PATHWAY_GENES <- gethsagenes(x=kegg_TARGET_PATHWAY)
names(kegg_TARGET_PATHWAY_GENES) <- kegg_TARGET_PATHWAY




