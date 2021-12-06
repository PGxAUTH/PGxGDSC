library(tidyverse)

Gain = 5.0
Loss = 0.5

GENE_NAMES <- rownames(CCLE_CN)

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

Cell_Lines <- colnames(CCLE_CN)

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

R_CN <- CCLE_CN %>% dplyr::select(dput(R_cell_lines))
S_CN <- CCLE_CN %>% dplyr::select(dput(S_cell_lines))

# Keep only GAIN & LOSS
R_CN_Gain <- R_CN %>% filter_all(all_vars(. > Gain))
S_CN_Gain <- S_CN %>% filter_all(all_vars(. > Gain))
R_GainCNgenes <- as.data.frame(rownames(R_CN_Gain))
colnames(R_GainCNgenes) <- "Genes"
S_GainCNgenes <- as.data.frame(rownames(S_CN_Gain))
colnames(S_GainCNgenes) <- "Genes"
Joined_GainCNgenes <- inner_join(R_GainCNgenes, S_GainCNgenes)[,1]

R_CN_Gain <- R_CN_Gain %>% mutate(Joined = rownames(R_CN_Gain) %in% Joined_GainCNgenes)
R_CN_Gain <- R_CN_Gain %>% filter(Joined == FALSE)
R_CN_Gain <- R_CN_Gain %>% dplyr::select(-Joined)
S_CN_Gain <- S_CN_Gain %>% mutate(Joined = rownames(S_CN_Gain) %in% Joined_GainCNgenes)
S_CN_Gain <- S_CN_Gain %>% filter(Joined == FALSE)
S_CN_Gain <- S_CN_Gain %>% dplyr::select(-Joined)

R_CN_Loss <- R_CN %>% filter_all(all_vars(. < Loss))
S_CN_Loss <- S_CN %>% filter_all(all_vars(. < Loss))
R_LossCNgenes <- as.data.frame(rownames(R_CN_Loss))
colnames(R_LossCNgenes) <- "Genes"
S_LossCNgenes <- as.data.frame(rownames(S_CN_Loss))
colnames(S_LossCNgenes) <- "Genes"
Joined_LossCNgenes <- inner_join(R_LossCNgenes, S_LossCNgenes)[,1]

R_CN_Loss <- R_CN_Loss %>% mutate(Joined = rownames(R_CN_Loss) %in% Joined_LossCNgenes)
R_CN_Loss <- R_CN_Loss %>% filter(Joined == FALSE)
R_CN_Loss <- R_CN_Loss %>% dplyr::select(-Joined)
S_CN_Loss <- S_CN_Loss %>% mutate(Joined = rownames(S_CN_Loss) %in% Joined_LossCNgenes)
S_CN_Loss <- S_CN_Loss %>% filter(Joined == FALSE)
S_CN_Loss <- S_CN_Loss %>% dplyr::select(-Joined)

R_CN_Gain_genes <- rownames(R_CN_Gain)
R_CN_Loss_genes <- rownames(R_CN_Loss)

R_CN_GainLoss <- rbind(R_CN_Gain, R_CN_Loss)
R_CN_GainLoss_genes <- rownames(R_CN_GainLoss)

S_CN_Gain_genes <- rownames(S_CN_Gain)
S_CN_Loss_genes <- rownames(S_CN_Loss)

S_CN_GainLoss <- rbind(S_CN_Gain, S_CN_Loss)
S_CN_GainLoss_genes <- rownames(S_CN_GainLoss)

# AVERAGE SCORE & SD
library(Matrix)
R_CN_AVERAGE <- rowMeans(R_CN, na.rm = TRUE)
R_CN_SD <- apply(R_CN[ ,1:ncol(R_CN)],1,sd)
R_CN <- R_CN %>% mutate(AVERAGE = R_CN_AVERAGE, SD = R_CN_SD)
S_CN_AVERAGE <- rowMeans(S_CN, na.rm = TRUE)
S_CN_SD <- apply(S_CN[ ,1:ncol(S_CN)],1,sd)
S_CN <- S_CN %>% mutate(AVERAGE = S_CN_AVERAGE, SD = S_CN_SD)

R_CN <- R_CN %>% mutate(GENE_NAME=GENE_NAMES)
S_CN <- S_CN %>% mutate(GENE_NAME=GENE_NAMES)

R_CN_GainLoss <- R_CN %>% filter(GENE_NAME %in% R_CN_GainLoss_genes == TRUE)
S_CN_GainLoss <- S_CN %>% filter(GENE_NAME %in% S_CN_GainLoss_genes == TRUE)

# Results
R_CN_GainLoss_left_joined <- left_join(R_CN_GainLoss, S_CN, by="GENE_NAME")
R_CN_GainLoss_left_joined <- R_CN_GainLoss_left_joined %>% transmute(GENE_NAME, AbsDiff=abs(AVERAGE.x - AVERAGE.y)) 

S_CN_GainLoss_left_joined <- left_join(S_CN_GainLoss, R_CN, by="GENE_NAME")
S_CN_GainLoss_left_joined <- S_CN_GainLoss_left_joined %>% transmute(GENE_NAME, AbsDiff=abs(AVERAGE.x - AVERAGE.y)) 


################################
R_CN_Gain_genes <- as.data.frame(str_split_fixed(R_CN_Gain_genes, "\\.\\.", 2))$V1
R_CN_Gain_genes <- as.data.frame(R_CN_Gain_genes)
colnames(R_CN_Gain_genes) <- "Genes"
R_CN_Gain_genes <- R_CN_Gain_genes %>% mutate("Cell Lines Class"="Resistant", "Source"="CCLE", "Type of analysis"="CNV", "Type of alteration"="Gain")

R_CN_Loss_genes <- as.data.frame(str_split_fixed(R_CN_Loss_genes, "\\.\\.", 2))$V1
R_CN_Loss_genes <- as.data.frame(R_CN_Loss_genes)
colnames(R_CN_Loss_genes) <- "Genes"
R_CN_Loss_genes <- R_CN_Loss_genes %>% mutate("Cell Lines Class"="Resistant", "Source"="CCLE", "Type of analysis"="CNV", "Type of alteration"="Loss")

R_CN_GainLoss_genes <- rbind(R_CN_Gain_genes, R_CN_Loss_genes)


S_CN_Gain_genes <- as.data.frame(str_split_fixed(S_CN_Gain_genes, "\\.\\.", 2))$V1
S_CN_Gain_genes <- as.data.frame(S_CN_Gain_genes)
colnames(S_CN_Gain_genes) <- "Genes"
S_CN_Gain_genes <- S_CN_Gain_genes %>% mutate("Cell Lines Class"="Sensitive", "Source"="CCLE", "Type of analysis"="CNV", "Type of alteration"="Gain")

S_CN_Loss_genes <- as.data.frame(str_split_fixed(S_CN_Loss_genes, "\\.\\.", 2))$V1
S_CN_Loss_genes <- as.data.frame(S_CN_Loss_genes)
colnames(S_CN_Loss_genes) <- "Genes"
S_CN_Loss_genes <- S_CN_Loss_genes %>% mutate("Cell Lines Class"="Sensitive", "Source"="CCLE", "Type of analysis"="CNV", "Type of alteration"="Loss")

S_CN_GainLoss_genes <- rbind(S_CN_Gain_genes, S_CN_Loss_genes)


CN_CCLE_RESULTS <- rbind(R_CN_GainLoss_genes, S_CN_GainLoss_genes)

CN_CCLE_RESULTS <- CN_CCLE_RESULTS %>% mutate_all(as.character)

core_geneset$CN_CCLE_RESULTS <- CN_CCLE_RESULTS




