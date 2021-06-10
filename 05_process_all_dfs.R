#This script does everything 04 does with 1 dataframe, but does it for all df in l_expr_tables

# l_expr_tables is opened
# remove gene id versions from dfs
# merge test data frame with pc table. Name it merg_expr_tables
# processes dfs in merge_expr_tables. This processing will create a new list containing dfs with the processed tables. It will be pro_expr_tables
# Get some visualization about the filtering is generated.
# Get some summary statistics about our data frame

# Returns:
#merg_expr_tables. The expression table test merged onto pc table
#pro_expr_tables. The expression table containing only unique expressed genes
#summary_visualization. visualization plot 


source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(skimr)

#-Open pc table
tryCatch(
  {
    pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
  },
  error = function(cond) {
    message(cond)
    message("Couldn't open the pc table. It should be in Data/")
  }
)

#---Open l_expr_tables
tryCatch(
  {
    l_expr_tables <- readRDS(file = "Data/l_expr_tables")
  },
  error = function(cond) {
    message(cond)
    message("Couldn't find l_expr_tables rds object. Should have been generated in 01")
  }
)

#---Remove gene id versions

l_expr_tables <- lapply(l_expr_tables, remove_gene_id_vers)

#---Merge experssion tables onto pc table

merg_expr_tables <- lapply(l_expr_tables, merge_tables, pc_table=pc)




