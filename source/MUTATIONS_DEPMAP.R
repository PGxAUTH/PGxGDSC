library(tidyverse)

# Select R & S cell lines
dataset <- dataset %>% mutate(IC50=exp(LN_IC50))
cell_lines <- dataset %>% filter(DRUG_NAME==drug_name, TCGA_DESC==tissue)
MaxConc <- as.numeric(cell_lines[1,"MAX_CONC"])
R <- cell_lines %>% filter(IC50 > MaxConc)
R <- R %>% filter(IC50>=unname(quantile(R$IC50, .75)))
R_cell_lines <- R$SANGER_MODEL_ID
S <- cell_lines %>% filter(IC50 <= MaxConc)
S <- S %>% filter(IC50<=unname(quantile(S$IC50, .25)))
S_cell_lines <- S$SANGER_MODEL_ID

# CHECK
Cell_Lines <- CCLE_mutations$Sanger_Model_ID

MatchR <- R_cell_lines %in% Cell_Lines 
MatchR <- data.frame(R_cell_lines, MatchR)
MatchR <- MatchR %>% filter(MatchR==TRUE)
R_cell_lines <- MatchR$R_cell_lines

MatchS <- S_cell_lines %in% Cell_Lines 
MatchS <- data.frame(S_cell_lines, MatchS)
MatchS <- MatchS %>% filter(MatchS==TRUE)
S_cell_lines <- MatchS$S_cell_lines

# STOP OR MOVE ON
{
  if (length(R_cell_lines)<=2 | length(S_cell_lines)<=2) {stop("There are less than three resistant or sensitive cell lines, the execution is stopped.")}
  
}

getlist <- function(cell_lines) {
  listofcell_lines <- list()
  for(i in 1:length(cell_lines)){
    df <- CCLE_mutations %>% filter(Sanger_Model_ID == cell_lines[i])
    listofcell_lines[[i]] <- df
    
  }
  return(listofcell_lines)
}

R_list <- getlist(cell_lines = R_cell_lines)
S_list <- getlist(cell_lines = S_cell_lines)

# R_common_mutations
R_common_mutations <- map(R_list, ~.$Hugo_Symbol) %>% 
  purrr::reduce(intersect)  
R_common_mutations_df <- as.data.frame(R_common_mutations)

# Check R_common_mutations %in% S
check1 <- function(S) {
  listofS <- list()
  for(i in 1:length(S)){
    vector <- R_common_mutations %in% S_list[[i]]$Hugo_Symbol
    listofS[[i]] <- vector
    
  }
  return(listofS)
}
S_check_list <- check1(S = S_list)

# Find R mutually exclusive mutated genes
S_check_df <- as.data.frame(S_check_list[])
S_check_df <- S_check_df %>% mutate(allFalse = (rowSums(S_check_df == TRUE) == 0 | rowSums(S_check_df == FALSE) == ncol(S_check_df)))

R_mutually_exclusive_mutated_genes <- cbind(R_common_mutations_df, S_check_df)
R_mutually_exclusive_mutated_genes <- R_mutually_exclusive_mutated_genes %>% filter(allFalse==TRUE)
R_mutually_exclusive_mutated_genes <- R_mutually_exclusive_mutated_genes[,1]
R_mutually_exclusive_mutated_genes

# S_common_mutations
S_common_mutations <- map(S_list, ~.$Hugo_Symbol) %>% 
  purrr::reduce(intersect)  
S_common_mutations_df <- as.data.frame(S_common_mutations)

# Check S_common_mutations %in% R
check2 <- function(R) {
  listofR <- list()
  for(i in 1:length(R)){
    vector <- S_common_mutations %in% R_list[[i]]$Hugo_Symbol
    listofR[[i]] <- vector
    
  }
  return(listofR)
}
R_check_list <- check2(R = R_list)

# Find S mutually exclusive mutated genes
R_check_df <- as.data.frame(R_check_list[])
R_check_df <- R_check_df %>% mutate(allFalse = (rowSums(R_check_df == TRUE) == 0 | rowSums(R_check_df == FALSE) == ncol(R_check_df)))

S_mutually_exclusive_mutated_genes <- cbind(S_common_mutations_df, R_check_df)
S_mutually_exclusive_mutated_genes <- S_mutually_exclusive_mutated_genes %>% filter(allFalse==TRUE)
S_mutually_exclusive_mutated_genes <- S_mutually_exclusive_mutated_genes[,1]
S_mutually_exclusive_mutated_genes

# PATHOGENIC Mutually Exclusive Mutated Genes (MEMG)
{
  if (length(R_mutually_exclusive_mutated_genes)>=1){
r_memg <- function(mutually_exclusive_mutated_genes) {
  listofmutually_exclusive_mutated_genes <- list()
  for(i in 1:length(mutually_exclusive_mutated_genes)){
    r_memg_list <- map(R_list, ~filter(.x, Hugo_Symbol == mutually_exclusive_mutated_genes[i]))
    listofmutually_exclusive_mutated_genes[[i]] <- r_memg_list
    
  }
  return(listofmutually_exclusive_mutated_genes)
}
R_MEMG_list <- r_memg(mutually_exclusive_mutated_genes = R_mutually_exclusive_mutated_genes)
  } else {
    R_MEMG_list <- list()
  }
}
  
{
  if (length(S_mutually_exclusive_mutated_genes)>=1){
s_memg <- function(mutually_exclusive_mutated_genes) {
  listofmutually_exclusive_mutated_genes <- list()
  for(i in 1:length(mutually_exclusive_mutated_genes)){
    s_memg_list <- map(S_list, ~filter(.x, Hugo_Symbol == mutually_exclusive_mutated_genes[i]))
    listofmutually_exclusive_mutated_genes[[i]] <- s_memg_list
    
  }
  return(listofmutually_exclusive_mutated_genes)
}
S_MEMG_list <- s_memg(mutually_exclusive_mutated_genes = S_mutually_exclusive_mutated_genes)
  } else {
    S_MEMG_list <- list()
  }
}

find_PATHOGENIC <- function(MEMG_list) {
  if (length(MEMG_list)>=1){
  vectorofMEMG_list <- as.character()
  for(i in 1:length(MEMG_list)){
    result <- as.logical(as.character(lapply(MEMG_list[[i]], function(x) if ("PATHOGENIC" %in% x$Variant_annotation==TRUE) TRUE else FALSE)))
    vectorofMEMG_list[[i]] <- ifelse(FALSE %in% result, "OTHER", MEMG_list[[i]][[1]]$Hugo_Symbol)
    
  }
  return(vectorofMEMG_list)
  } else {
    return(as.character(vector()))
  }
}

R_MEMG_PATHOGENIC <- find_PATHOGENIC(MEMG_list = R_MEMG_list)
R_MEMG_PATHOGENIC = R_MEMG_PATHOGENIC[R_MEMG_PATHOGENIC!= "OTHER"]
R_MEMG_PATHOGENIC

S_MEMG_PATHOGENIC <- find_PATHOGENIC(MEMG_list = S_MEMG_list)
S_MEMG_PATHOGENIC = S_MEMG_PATHOGENIC[S_MEMG_PATHOGENIC!= "OTHER"]
S_MEMG_PATHOGENIC



################################
R_mutually_exclusive_mutated_genes <- as.data.frame(R_mutually_exclusive_mutated_genes)
colnames(R_mutually_exclusive_mutated_genes) <- "Genes"
R_mutually_exclusive_mutated_genes <- R_mutually_exclusive_mutated_genes %>% mutate("Cell Lines Class"="Resistant", "Source"="CCLE", "Type of analysis"="Mutation", "Type of alteration"=ifelse(`Genes` %in% R_MEMG_PATHOGENIC==TRUE, "Damaging", "Other"))

S_mutually_exclusive_mutated_genes <- as.data.frame(S_mutually_exclusive_mutated_genes)
colnames(S_mutually_exclusive_mutated_genes) <- "Genes"
S_mutually_exclusive_mutated_genes <- S_mutually_exclusive_mutated_genes %>% mutate("Cell Lines Class"="Sensitive", "Source"="CCLE", "Type of analysis"="Mutation", "Type of alteration"=ifelse(`Genes` %in% S_MEMG_PATHOGENIC==TRUE, "Damaging", "Other"))

MUTATIONS_CCLE_RESULTS <- rbind(R_mutually_exclusive_mutated_genes, S_mutually_exclusive_mutated_genes)

MUTATIONS_CCLE_RESULTS <- MUTATIONS_CCLE_RESULTS %>% mutate_all(as.character)

core_geneset$MUTATIONS_CCLE_RESULTS <- MUTATIONS_CCLE_RESULTS


