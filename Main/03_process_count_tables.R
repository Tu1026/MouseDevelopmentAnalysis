#This script does the following
#1) Open all of the count tables in Data/count_tables as a list named l_expr_tables
#2) Processes a single count table
#3) Generate Count-sample matrix

source('functions.R')

library(tidyverse)
library(assertthat)
library(stringr)
library(skimr)
library(Hmisc)
library(corrplot)
library(dplyr)


#---------------------------------------------------------------------------
########################## 1) Opening Count Tables
#--------------------------------------------------------------------------

# Opens all of the count tables in Data/count_tables as a list named l_expr_tables

#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")

#---------------------------------------------------------------------------
# Global variables
#--------------------------------------------------------------------------
target_directory <- "Data/Count_tables"

#---------------------------------------------------------------------------
# Setting the order that the count tables are opened in
#--------------------------------------------------------------------------
# We are opening the count tables into a list. We want this list to be
# ordered in the same order as the meta_data table. This is to 
# represent the samples in order of develempental times

ordered_meta_data_id_names <- paste0(meta_data$id,".tsv")

#---------------------------------------------------------------------------
# Open the count tables
#--------------------------------------------------------------------------
# Open count tables and save inside a list. List should be ordered
# In the order of the meta_data table.
# Remove the .tsv after the names once the files are read

l_expr_tables <- lapply(ordered_meta_data_id_names,
                        open_expr_table,
                        target_directory)

names(l_expr_tables) <- str_replace(ordered_meta_data_id_names, ".tsv","")



#---------------------------------------------------------------------------
########################## 2) Generate  sample - expression data matrix
#--------------------------------------------------------------------------

# This script processes ALL count tables within l_expr_tables. 
# remove gene id versions from test data frames
# merge test data frame with pc table.
# Data frames are used to generate a  sample - expression data matrix
# sample - expression data matrix is saved


#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

# l_expr_tables is already opened: It is a list containing all of the count tables
# pc_sub is already generated: It is a processed pc table

#---Open meta_data
meta_data <- read.delim("Data/complete_meta_data.tsv",
                        stringsAsFactors = FALSE,
                        sep = "\t")

#---Open PC Table
pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv",
                 stringsAsFactors = FALSE)

#---Remove nondistinct Symbols and Gene_IDs from pc table
# and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  dplyr::select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)

#---------------------------------------------------------------------------
# Prep PC and l_expr_tables for merging with each other, then merge
#--------------------------------------------------------------------------

#---Run remove_gene_id_vers on count tables to remove gene id versions
stopifnot("Gene_ID" %in% colnames(l_expr_tables[[1]]))
l_expr_tables_noid <- lapply(l_expr_tables, remove_gene_id_vers)

#---merge pc_sub table with all count tables
l_expr_merged <- lapply(l_expr_tables_noid, merge_count_pc, pc_sub)

#---------------------------------------------------------------------------
# Get only Protein Coding genes from merged expression - pc data frames. 
#---------------------------------------------------------------------------
# Now that our count tables have Protein Symbols, we can 
# Use the presence of a symbol to filter for protein coding genes 
l_expr_symbols <- lapply(l_expr_merged, get_genes_with_symbols)

#---------------------------------------------------------------------------
# Creating the sample - expression data matrix
#--------------------------------------------------------------------------
#---Test to see that each frame in l_expr_symbols has the same order of symbols
# If this is true, then we can just iteratively add the lists onto a matrix

stopifnot(test_orders(l_expr_symbols))


#---Create empty expression matrix and fill with expression data
count_matrix <- matrix(data = 0,
                            nrow = nrow(l_expr_symbols$ENCFF227HKF),
                            ncol = length(l_expr_symbols))

rownames(count_matrix) <- l_expr_symbols$ENCFF227HKF$Symbol
colnames(count_matrix) <- names(l_expr_symbols)

for (colname in colnames(count_matrix)) { 
  count_matrix[ , colname] <- l_expr_symbols[[colname]][["TPM"]]
}

#---Save expression matrix as rds

saveRDS(object = count_matrix, file = "Data/pc_count_matrix.rds")

#---------------------------------------------------------------------------
# Building an AVG Matrix
#---------------------------------------------------------------------------
# As seen in meta_data, every mouse developmental stage has 2 samples assosiated
# with it. We want to create a new matrix with those replicates averaged so
# we can easily perform analysis on how our counts change over the dev stages
# Build a new matrix containing the averages of counts across replicates
# Called avg_count_matrix

uniq_dev_stages <- unique(meta_data$dev_stage)

avg_count_matrix <- matrix(data = 0,
                           nrow = nrow(count_matrix),
                           ncol = ncol(count_matrix)/2)
rownames(avg_count_matrix) <- rownames(count_matrix)
colnames(avg_count_matrix) <- uniq_dev_stages 

for (mydev_stage in uniq_dev_stages) { 
  meta_subset <- filter(meta_data, meta_data$dev_stage == mydev_stage)
  avg <- rowMeans(count_matrix[,meta_subset$id])
  avg_count_matrix[,mydev_stage] <- avg
}

#---Check to confirm the averaging functioned properly

for (check_number in 1:5) {
  check_name <- sample(rownames(avg_count_matrix), size = 1)
  
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,1:2]),
                                  avg_count_matrix[check_name, 1]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,3:4]),
                                  avg_count_matrix[check_name, 2]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,9:10]),
                                  avg_count_matrix[check_name, 5]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,15:16]),
                                  avg_count_matrix[check_name, 8]))
}

#---Save avg expression matrix as rds

saveRDS(object = avg_count_matrix, file = "Data/avg_pc_count_matrix.rds")


#---------------------------------------------------------------------------
# Creating a difference count matrix
#---------------------------------------------------------------------------
# We want to identify which genes in a pair of replicates differ the most
# in count value.


uniq_dev_stages <- unique(meta_data$dev_stage)

diff_count_matrix <- matrix(0, 
                            nrow = nrow(avg_count_matrix),
                            ncol = ncol(avg_count_matrix))
rownames(diff_count_matrix) <- rownames(avg_count_matrix)
colnames(diff_count_matrix) <- colnames(avg_count_matrix)

for (mydev_stage in uniq_dev_stages) { 
  meta_subset <- filter(meta_data, meta_data$dev_stage == mydev_stage)
  stopifnot(length(meta_subset$id) == 2)
  diff_count_matrix[, mydev_stage] <-  abs(count_matrix[ , meta_subset$id[1]] - count_matrix[ , meta_subset$id[2]])
}


#---Check to confirm the averaging functioned properly

for (check_number in 1:5) {
  check_name <- sample(rownames(diff_count_matrix), size = 1)
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 1] - count_matrix[check_name, 2])
                                  , diff_count_matrix[check_name, 1]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 3] - count_matrix[check_name, 4])
                                  , diff_count_matrix[check_name, 2]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 9] - count_matrix[check_name, 10])
                                  , diff_count_matrix[check_name, 5]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 15] - count_matrix[check_name, 16])
                                  , diff_count_matrix[check_name, 8]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 13] - count_matrix[check_name, 14])
                                  , diff_count_matrix[check_name, 7]))
  
}


#---Save avg expression matrix as rds

saveRDS(object = diff_count_matrix, file = "Data/diff_pc_count_matrix.rds")



