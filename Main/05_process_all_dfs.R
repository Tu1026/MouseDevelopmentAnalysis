#This script processes ALL count tables within l_expr_tables. 
# l_expr_tables is opened,
# remove gene id versions from test data frames
# merge test data frame with pc table.
# Data frames are processed. Generating l_with_symbols. Which contains the protein coding genes for each count table

source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(skimr)

#---Open l_expr_tables

if (file.exists(file = "Data/l_expr_tables.rds")) {
  message("Found l_expr_tables.rds, opening...")
  l_expr_tables <- readRDS(file = "Data/l_expr_tables.rds")
} else {
  message("Couldn't find l_expr_tables rds object. Should have been generated in 01")
}



#---Open pc table

if (file.exists(file = "Data/ensembl_mouse_protein_coding_104.tsv")) {
  message("protein coding table found, opening...")
  pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
} else {
  message("Couldn't find pc table, it should be been downloaded in ")
}

#---Remove nondistinct from pc table and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  select(Gene_ID, Symbol)

#---Runn remove_gene_id_vers on count tables to remove idversions

l_expr_tables_noid <- lapply(l_expr_tables, remove_gene_id_vers)

#---merge pc table with all count tables

merge_count_pc <- function(count_table, pc_sub) { 
  #merge a count table onto the pc table
  merged_table <- count_table %>%
    dplyr::rename(Gene_ID = gene_id) %>% 
    left_join(y = pc_sub, by = "Gene_ID") %>% 
    relocate(Symbol, .after = Gene_ID)
  return (merged_table)
  }

l_expr_merged <- lapply(l_expr_tables_noid, merge_count_pc, pc_sub)

#---------------------------------------------------------------------------
# This part of the  script processes the merged expression - pc data frames. 
#---------------------------------------------------------------------------
process_merged_table <- function(merged_table) {
      #Do a bunch of processes, returning with_symbols, the data frame that has protein coding counts with symbols
      
    #---get a data frame of all the empty symbols
    
    empty_symbols <- merged_table%>%
      filter(is.na(Symbol))
    
    #---Filter out transcripts that do have a protein symbol
    
    with_symbols <- merged_table %>%
      filter(!(is.na(Symbol)))
    
    #---Check to make sure empty_symbols + with_symbols is equal to the total amount of rows in test_merged
    stopifnot(nrow(merged_table) == nrow(empty_symbols)+nrow(with_symbols))
    
    #---Check that all of your transcriopts in with_symbols do indeed have symbols
    stopifnot(!(all(is.na(with_symbols$Symbol))))
    
    #---Check that all of your transcriopts in empty_symbols do NOT have symbols
    stopifnot((all(is.na(empty_symbols$Symbol))))
    
    #---Make sure the number of distinct genes was equal to the num of groups before
    n_groups <- merged_table%>%
      filter(!is.na(Symbol))%>%
      group_by(Gene_ID)%>%
      tally()%>%
      nrow()
    stopifnot(nrow(with_symbols) == n_groups)
    
    return (with_symbols)
}

l_with_symbols <- lapply(l_expr_merged, process_merged_table)

# test_frame[duplicated(test_frame$gene_id),]
# test_merged[duplicated(test_merged$Symbol),]
# filter(test_merged, test_merged$Symbol=="Pcdha11")
# test_with_symbols <- test_merged %>%
#   filter(!(is.na(Symbol)))
# glimpse(test_with_symbols[duplicated(test_with_symbols$Symbol),])
# glimpse(pc_sub[duplicated(pc_sub$Symbol),])

#---Test to see that each frame in l_with_symbols has the same order of symbolss
test_orders <- function(l_with_symbols) { 
  reference_frame <- l_with_symbols$ENCFF227HKF$Symbol 
  for (i in 1:length(l_with_symbols)) {
    test_frame <- l_with_symbols[[i]]$Symbol
    stopifnot(all(reference_frame == test_frame))
  }
  message("The order of all symbols within l_with_symbols is the same")
  return(TRUE)
}


create_data_matrix <- function(l_with_symbols, dimnames) {
  matrix = (l_with_symbols[[1]]$TPM)
            
  for (i in 2:length(l_with_symbols)) {
    matrix <- cbind(matrix, l_with_symbols[[i]]$TPM)
  } 
  dimnames(matrix) <- dimnames
  return (matrix)
}

col_names <- names(l_with_symbols)
row_names <- l_with_symbols[[1]]$Symbol
my_dimnames <- list(row_names, col_names) 

if (test_orders(l_with_symbols)) {
  expression_matrix <- create_data_matrix(l_with_symbols, my_dimnames)
}

#---Save expression matrix as tab delimited

saveRDS(object = expression_matrix, file = "Data/protein_expression_matrix")

