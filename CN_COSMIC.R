library(tidyverse)

# Select R & S cell lines
dataset <- dataset %>% mutate(IC50=exp(LN_IC50))
cell_lines <- dataset %>% filter(DRUG_NAME==drug_name, TCGA_DESC==tissue)
MaxConc <- as.numeric(cell_lines[1,"MAX_CONC"])
R <- cell_lines %>% filter(IC50 > MaxConc)
R <- R %>% filter(IC50>=unname(quantile(R$IC50, .75)))
R_cell_lines <- R$CELL_LINE_NAME
S <- cell_lines %>% filter(IC50 <= MaxConc)
S <- S %>% filter(IC50<=unname(quantile(S$IC50, .25)))
S_cell_lines <- S$CELL_LINE_NAME

# CHECK
Cell_Lines <- COSMIC_CNVs$SAMPLE_NAME

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
    df <- COSMIC_CNVs %>% filter(SAMPLE_NAME == cell_lines[i])
    listofcell_lines[[i]] <- df
    
  }
  return(listofcell_lines)
}

R_list <- getlist(cell_lines = R_cell_lines)
S_list <- getlist(cell_lines = S_cell_lines)

# R_common_CNVs
R_common_CNVs <- map(R_list, ~.$Gene_MUT_TYPE) %>% 
  purrr::reduce(intersect)  
R_common_CNVs_df <- as.data.frame(R_common_CNVs)

# Check R_common_CNVs %in% S
check1 <- function(S) {
  listofS <- list()
  for(i in 1:length(S)){
    vector <- R_common_CNVs %in% S_list[[i]]$Gene_MUT_TYPE
    listofS[[i]] <- vector
    
  }
  return(listofS)
}
S_check_list <- check1(S = S_list)

# Find R mutually exclusive CNVs
S_check_df <- as.data.frame(S_check_list[])
S_check_df <- S_check_df %>% mutate(allFalse = (rowSums(S_check_df == TRUE) == 0 | rowSums(S_check_df == FALSE) == ncol(S_check_df)))

R_mutually_exclusive_CNVs <- cbind(R_common_CNVs_df, S_check_df)
R_mutually_exclusive_CNVs <- R_mutually_exclusive_CNVs %>% filter(allFalse==TRUE)
R_mutually_exclusive_CNVs <- R_mutually_exclusive_CNVs[,1]
R_mutually_exclusive_CNVs

# S_common_CNVs
S_common_CNVs <- map(S_list, ~.$Gene_MUT_TYPE) %>% 
  purrr::reduce(intersect)  
S_common_CNVs_df <- as.data.frame(S_common_CNVs)

# Check S_common_mutations %in% R
check2 <- function(R) {
  listofR <- list()
  for(i in 1:length(R)){
    vector <- S_common_CNVs %in% R_list[[i]]$Gene_MUT_TYPE
    listofR[[i]] <- vector
    
  }
  return(listofR)
}
R_check_list <- check2(R = R_list)

# Find S mutually exclusive CNVs
R_check_df <- as.data.frame(R_check_list[])
R_check_df <- R_check_df %>% mutate(allFalse = (rowSums(R_check_df == TRUE) == 0 | rowSums(R_check_df == FALSE) == ncol(R_check_df)))

S_mutually_exclusive_CNVs <- cbind(S_common_CNVs_df, R_check_df)
S_mutually_exclusive_CNVs <- S_mutually_exclusive_CNVs %>% filter(allFalse==TRUE)
S_mutually_exclusive_CNVs <- S_mutually_exclusive_CNVs[,1]
S_mutually_exclusive_CNVs


################################
R_mutually_exclusive_CNVs <- as.data.frame(R_mutually_exclusive_CNVs)
colnames(R_mutually_exclusive_CNVs) <- "Genes"
R_mutually_exclusive_CNVs <- R_mutually_exclusive_CNVs %>% mutate("Genes"=as.data.frame(str_split_fixed(R_mutually_exclusive_CNVs$`Genes`, ", ", 2))$V1, "Cell Lines Class"="Resistant", "Source"="COSMIC", "Type of analysis"="CNV", "Type of alteration"=as.data.frame(str_split_fixed(R_mutually_exclusive_CNVs$`Genes`, ", ", 2))$V2)

S_mutually_exclusive_CNVs <- as.data.frame(S_mutually_exclusive_CNVs)
colnames(S_mutually_exclusive_CNVs) <- "Genes"
S_mutually_exclusive_CNVs <- S_mutually_exclusive_CNVs %>% mutate("Genes"=as.data.frame(str_split_fixed(S_mutually_exclusive_CNVs$`Genes`, ", ", 2))$V1, "Cell Lines Class"="Sensitive", "Source"="COSMIC", "Type of analysis"="CNV", "Type of alteration"=as.data.frame(str_split_fixed(S_mutually_exclusive_CNVs$`Genes`, ", ", 2))$V2)

CN_COSMIC_RESULTS <- rbind(R_mutually_exclusive_CNVs, S_mutually_exclusive_CNVs)

CN_COSMIC_RESULTS <- CN_COSMIC_RESULTS %>% mutate_all(as.character)

core_geneset$CN_COSMIC_RESULTS <- CN_COSMIC_RESULTS


