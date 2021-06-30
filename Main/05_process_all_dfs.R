#This script processes ALL count tables within l_expr_tables. 
# l_expr_tables is opened,
# remove gene id versions from test data frames
# merge test data frame with pc table.
# Data frames are processed. Generating l_expr_symbols. Which contains the protein coding genes for each count table

source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(skimr)
library(Hmisc)
library(corrplot)
#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

#open l_expr_tables.rds
#An RDS object containing all of the expression tables within a List 

l_expr_tables <- readRDS(file = "Data/l_expr_tables.rds")

#---Open pc table

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)

#---Open meta_data

meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")


#---------------------------------------------------------------------------
# Prep PC and l_expr_tables for merging with eachother, then merge
#--------------------------------------------------------------------------

#---Remove nondistinct from pc table and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)

#---Run remove_gene_id_vers on count tables to remove gene id versions

l_expr_tables_noid <- lapply(l_expr_tables, remove_gene_id_vers)

#---merge pc table with all count tables

merge_count_pc <- function(count_table, pc_sub) { 
  #merge a count table onto the pc table
  merged_table <- count_table %>%
    left_join(y = pc_sub, by = "Gene_ID") %>% 
    relocate(Symbol, .after = Gene_ID)
  return (merged_table)
  }

l_expr_merged <- lapply(l_expr_tables_noid, merge_count_pc, pc_sub)

#---------------------------------------------------------------------------
# Get only Protein Coding genes from merged expression - pc data frames. 
#---------------------------------------------------------------------------
get_genes_with_symbols <- function(merged_table) {
    #---Return a table containing the transcripts that have a protein symbol
    
    with_symbols <- merged_table %>%
      filter(!(is.na(Symbol)))
    
    return (with_symbols)
}

l_expr_symbols <- lapply(l_expr_merged, get_genes_with_symbols)

#---------------------------------------------------------------------------
# Creating the sample - expression data matrix
#---------------------------------------------------------------------------

#---Test to see that each frame in l_expr_symbols has the same order of symbolss
test_orders <- function(l_expr_symbols) { 
  reference_frame <- l_expr_symbols$ENCFF227HKF$Symbol 
  for (i in 1:length(l_expr_symbols)) {
    test_frame <- l_expr_symbols[[i]]$Symbol
    stopifnot(all(reference_frame == test_frame))
  }
  message("The order of all symbols within l_expr_symbols is the same")
  return(TRUE)
}

stopifnot(test_orders(l_expr_symbols))

#---Create empty expression matrix and fill with expression data
expression_matrix <- matrix(data = 0,
                            nrow = nrow(l_expr_symbols$ENCFF227HKF),
                            ncol = length(l_expr_symbols))
rownames(expression_matrix) <- l_expr_symbols$ENCFF227HKF$Symbol
colnames(expression_matrix) <- names(l_expr_symbols)

for (i in 1:length(l_expr_symbols)) {
  expression_matrix[,i] <- l_expr_symbols[[i]]$TPM
}

#---Save expression matrix as rds

saveRDS(object = expression_matrix, file = "Data/protein_expression_matrix.rds")

#---------------------------------------------------------------------------
# Generating Sample Correlation Matrix
#---------------------------------------------------------------------------
# Generating correlation matrix to see correlation between sample's expression
# Also generate correlation matrix's p-values

correlation_matrix_rcorr <- rcorr(expression_matrix, type = "spearman")
correlation_matrix <- correlation_matrix_rcorr$r
correlation_matrix_pvalues <- correlation_matrix_rcorr$P

#---Generating heatmap for sample correlation matrix

heatmap <- heatmap(x = correlation_matrix, sym = TRUE, Rowv = NA, Colv = NA, revC= TRUE)

#---------------------------------------------------------------------------
# Compressing Metadata 
#---------------------------------------------------------------------------
# Not sure if there are betters ways to extract specific columns from matrices,
# But this is the best I could do
# I'm literally just trying to figure out how to link the replicates onto the expression matrix columns

colnames <- colnames(expression_matrix)


meta_data_group <- meta_data %>%
  select(id, dev_stage) %>%
  group_by(dev_stage)%>%
  dplyr::summarize(ids = str_split(paste(id, collapse = ","), pattern = ","))
rownames(meta_data_group) <- meta_data_group$dev_stage

indexes <- list()
for (i in meta_data_group$dev_stage) {
  row <- meta_data_group[i,]
  names <- unlist(row$ids)
  cur_indexes <- which(colnames == names)
  indexes <- append(indexes,list(cur_indexes))
}

indexes_almost_merged <- do.call(rbind, indexes)
meta_data_group <- meta_data_group %>%
  mutate(index1 = indexes_almost_merged[,1],
         index2 = indexes_almost_merged[,2])

#---------------------------------------------------------------------------
# Comparing reads between replicates
#---------------------------------------------------------------------------


#-------doing for 1
               
e105 <- meta_data_group["e10.5",]
e105 <- unlist(e105$ids)
indexes <- which(colnames == e105)

extract_matrix_columns <- function(index, expression_matrix) {
  return (expression_matrix[,index])
}

e105_expression_matrix <- do.call(cbind, lapply(indexes, 
                                                extract_matrix_columns,
                                                expression_matrix))

replicate_differences <- abs(e105_expression_matrix[,1] - e105_expression_matrix[,2])

#-------Doing for all count tables


master_diff_matrix <- matrix(data = 0,
                             nrow = nrow(expression_matrix), 
                             ncol = nrow(meta_data_group),   
                             dimnames =  list(rownames(expression_matrix), meta_data_group$dev_stage))
for (i in 1:nrow(meta_data_group)) {
  row <- meta_data_group[i,]
  matrixcol1 <- expression_matrix[,row$index1]
  matrixcol2 <- expression_matrix[,row$index2]
  replicate_diff <- abs(matrixcol1-matrixcol2)
  master_diff_matrix[,i] <-  replicate_diff
}


#---------------------------------------------------------------------------
# Averaging Replicates
#---------------------------------------------------------------------------
# again I'm fkinda floundering here. I bet there's better ways 

#----doing for all count tables

master_avg_matrix <- matrix(data = 0 , 
                            ncol = nrow(meta_data_group),
                            nrow = nrow(expression_matrix), 
                            dimnames =  list(rownames(expression_matrix), meta_data_group$dev_stage))

for (i in 1:nrow(meta_data_group)) {
  row <- meta_data_group[i,]
  matrixcol1 <- expression_matrix[,row$index1]
  matrixcol2 <- expression_matrix[,row$index2]
  matrix_combined <- cbind(matrixcol1, matrixcol2)
  replicate_avg <-  round(rowMeans(matrix_combined), digits = 3)
  master_avg_matrix[,i] <-  replicate_avg
}





