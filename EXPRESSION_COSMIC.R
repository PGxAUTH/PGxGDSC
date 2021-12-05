library(tidyverse)

Z_SCORE_over = 2.001
Z_SCORE_under = -2.001

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
Cell_Lines <- colnames(COSMIC_exp)

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
R_exp <- COSMIC_exp %>% dplyr::select(dput(R_cell_lines))
suppressWarnings (rownames(R_exp)  <-  as.character(COSMIC_exp$GENE_NAME))
S_exp <- COSMIC_exp %>% dplyr::select(dput(S_cell_lines))
suppressWarnings (rownames(S_exp)  <-  as.character(COSMIC_exp$GENE_NAME))

GENE_NAMES <- COSMIC_exp$GENE_NAME

# Keep only OVER & UNDER
R_exp_Over <- R_exp %>% filter_all(all_vars(. >= Z_SCORE_over))
S_exp_Over <- S_exp %>% filter_all(all_vars(. >= Z_SCORE_over))
R_overexpgenes <- as.data.frame(rownames(R_exp_Over))
colnames(R_overexpgenes) <- "Genes"
S_overexpgenes <- as.data.frame(rownames(S_exp_Over))
colnames(S_overexpgenes) <- "Genes"
Joined_overexpgenes <- inner_join(R_overexpgenes, S_overexpgenes)[,1]

R_exp_Over <- R_exp_Over %>% mutate(Joined = rownames(R_exp_Over) %in% Joined_overexpgenes)
R_exp_Over <- R_exp_Over %>% filter(Joined == FALSE)
R_exp_Over_genes <- rownames(R_exp_Over)
R_exp_Over <- R_exp_Over %>% dplyr::select(-Joined)
suppressWarnings (rownames(R_exp_Over) <- R_exp_Over_genes)
S_exp_Over <- S_exp_Over %>% mutate(Joined = rownames(S_exp_Over) %in% Joined_overexpgenes)
S_exp_Over <- S_exp_Over %>% filter(Joined == FALSE)
S_exp_Over_genes <- rownames(S_exp_Over)
S_exp_Over <- S_exp_Over %>% dplyr::select(-Joined)
suppressWarnings (rownames(S_exp_Over) <- S_exp_Over_genes)

R_exp_Under <- R_exp %>% filter_all(all_vars(. <= Z_SCORE_under))
S_exp_Under <- S_exp %>% filter_all(all_vars(. <= Z_SCORE_under))
R_underexpgenes <- as.data.frame(rownames(R_exp_Under))
colnames(R_underexpgenes) <- "Genes"
S_underexpgenes <- as.data.frame(rownames(S_exp_Under))
colnames(S_underexpgenes) <- "Genes"
Joined_underexpgenes <- inner_join(R_underexpgenes, S_underexpgenes)[,1]

R_exp_Under <- R_exp_Under %>% mutate(Joined = rownames(R_exp_Under) %in% Joined_underexpgenes)
R_exp_Under <- R_exp_Under %>% filter(Joined == FALSE)
R_exp_Under_genes <- rownames(R_exp_Under)
R_exp_Under <- R_exp_Under %>% dplyr::select(-Joined)
suppressWarnings (rownames(R_exp_Under) <- R_exp_Under_genes)
S_exp_Under <- S_exp_Under %>% mutate(Joined = rownames(S_exp_Under) %in% Joined_underexpgenes)
S_exp_Under <- S_exp_Under %>% filter(Joined == FALSE)
S_exp_Under_genes <- rownames(S_exp_Under)
S_exp_Under <- S_exp_Under %>% dplyr::select(-Joined)
suppressWarnings (rownames(S_exp_Under) <- S_exp_Under_genes)

R_exp_OverUnder <- rbind(R_exp_Over, R_exp_Under)
R_exp_OverUnder_GENE_NAMES <- rownames(R_exp_OverUnder)

S_exp_OverUnder <- rbind(S_exp_Over, S_exp_Under)
S_exp_OverUnder_GENE_NAMES <- rownames(S_exp_OverUnder)

# AVERAGE Z SCORE & SD
library(Matrix)
R_exp_AVERAGE <- rowMeans(R_exp, na.rm = TRUE)
R_exp_SD <- apply(R_exp[ ,1:ncol(R_exp)],1,sd)
R_exp <- R_exp %>% mutate(GENE_NAME = GENE_NAMES, AVERAGE_Z_SCORE = R_exp_AVERAGE, SD = R_exp_SD)
S_exp_AVERAGE <- rowMeans(S_exp, na.rm = TRUE)
S_exp_SD <- apply(S_exp[ ,1:ncol(S_exp)],1,sd)
S_exp <- S_exp %>% mutate(GENE_NAME = GENE_NAMES, AVERAGE_Z_SCORE = S_exp_AVERAGE, SD = S_exp_SD)

R_exp_OverUnder <- R_exp %>% filter(GENE_NAME %in% R_exp_OverUnder_GENE_NAMES == TRUE)
S_exp_OverUnder <- S_exp %>% filter(GENE_NAME %in% S_exp_OverUnder_GENE_NAMES == TRUE)

# Results
R_exp_OverUnder_left_joined <- left_join(R_exp_OverUnder, S_exp, by="GENE_NAME")
R_exp_OverUnder_left_joined <- R_exp_OverUnder_left_joined %>% transmute(GENE_NAME, AbsDiff=abs(AVERAGE_Z_SCORE.x - AVERAGE_Z_SCORE.y)) 

S_exp_OverUnder_left_joined <- left_join(S_exp_OverUnder, R_exp, by="GENE_NAME")
S_exp_OverUnder_left_joined <- S_exp_OverUnder_left_joined %>% transmute(GENE_NAME, AbsDiff=abs(AVERAGE_Z_SCORE.x - AVERAGE_Z_SCORE.y)) 

################################
R_exp_Over_genes <- as.data.frame(R_exp_Over_genes)
colnames(R_exp_Over_genes) <- "Genes"
R_exp_Over_genes <- R_exp_Over_genes %>% mutate("Cell Lines Class"="Resistant", "Source"="COSMIC", "Type of analysis"="Expression", "Type of alteration"="Over")

R_exp_Under_genes <- as.data.frame(R_exp_Under_genes)
colnames(R_exp_Under_genes) <- "Genes"
R_exp_Under_genes <- R_exp_Under_genes %>% mutate("Cell Lines Class"="Resistant", "Source"="COSMIC", "Type of analysis"="Expression", "Type of alteration"="Under")

R_exp_OverUnder_genes <- rbind(R_exp_Over_genes, R_exp_Under_genes)


S_exp_Over_genes <- as.data.frame(S_exp_Over_genes)
colnames(S_exp_Over_genes) <- "Genes"
S_exp_Over_genes <- S_exp_Over_genes %>% mutate("Cell Lines Class"="Sensitive", "Source"="COSMIC", "Type of analysis"="Expression", "Type of alteration"="Over")

S_exp_Under_genes <- as.data.frame(S_exp_Under_genes)
colnames(S_exp_Under_genes) <- "Genes"
S_exp_Under_genes <- S_exp_Under_genes %>% mutate("Cell Lines Class"="Sensitive", "Source"="COSMIC", "Type of analysis"="Expression", "Type of alteration"="Under")

S_exp_OverUnder_genes <- rbind(S_exp_Over_genes, S_exp_Under_genes)


EXPRESSION_COSMIC_RESULTS <- rbind(R_exp_OverUnder_genes, S_exp_OverUnder_genes)

EXPRESSION_COSMIC_RESULTS <- EXPRESSION_COSMIC_RESULTS %>% mutate_all(as.character)

core_geneset$EXPRESSION_COSMIC_RESULTS <- EXPRESSION_COSMIC_RESULTS

