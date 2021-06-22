#This script opens all of the count tables in Data/count_tables as a list.

#Returns:
#It then saves this list as an rds object named l_expr_tables. 


source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(stringr)

#---functions
#Open a list of tables using the return of create_list_of_tables

# open_expr_tables <- function(count_dir) {
#   #
#   #creates a list containg all expression table dataframes
#   #
#   #@param: String of the directory name where data frames live. This is passed into create_names_of_tables.
#   #example: count_dir <- "Data/Count_tables"
#   #@return: A list containing all of the expression table dataframes. They are named according to their name within their directory
#   names_of_tables <- create_names_of_tables(count_dir)
#   l_expr_tables <- lapply(names_of_tables,read.delim)
#   
#   names <- list.files(count_dir)
#   names <- str_replace(names, pattern = ".tsv", replacement = "")
#   names(l_expr_tables) <- names
#   
#   return (l_expr_tables)
# }


open_expr_table <- function(expression_table.tsv, target_directory, column_names = c("gene_id", "transcript_id.s." , "length", "effective_length","expected_count", "TPM", "FPKM",
                                                                                     "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "TPM_ci_lower_bound", "TPM_ci_upper_bound",                                                                                   "FPKM_ci_lower_bound", "FPKM_ci_upper_bound")) {
  #open a single count experession table. 
  #expression_table.tsv = the name of the expression table you wish to open as string
  #target_directory = the directory where the expression table is saved, as string
  #column_names = the column names you wish to open the expression table as. Setting this as a vector will ovverwite the default names, a vector strings
  
  tryCatch(
    {
    table <- read.delim(file = paste0(target_directory,"/",expression_table.tsv), sep = "\t", col.names = column_names)
    return(table)
    },
    error = function(cond){
      message(cond)
      return(NA)
    }
  )
}


#---Rename list intexes to table names
rename_count_tables <- function(l_expr_tables, list_of_expression_tables) {
  names <- str_replace(list_of_expression_tables, pattern = ".tsv", replacement = "")
  names(l_expr_tables) <- names
  return(l_expr_tables)
}

#---function to open all tables within 1 list
open_expr_tables <- function(target_directory, list_of_expression_tables) {
  #Open up all of the expression tables in list_of_expression_tables. All tables must be tsvs, and they must all be in target_directory
  #return l_expr_tables, which is also saved as an rds in target_directory
  
  l_expr_tables <- lapply(list_of_expression_tables, open_expr_table, target_directory = target_directory)
  l_expr_tables <- rename_count_tables(l_expr_tables, list_of_expression_tables)
  saveRDS(l_expr_tables, file = "Data/l_expr_tables.rds")
  return(l_expr_tables)
}

#---implementation
target_directory <- "Data/Count_tables"
list_of_expression_tables <- list.files(target_directory)
l_expr_tables <- open_expr_tables(target_directory, list_of_expression_tables = list_of_expression_tables)


# #---directory where the count tables are at
# count_dir <- "Data/Count_tables"
# 
# 
# if (!(file.exists("Data/l_expr_tables.rds"))) {
#   message ("l_expr_tables has not been created yet. Creating ... ")
#   l_expr_tables <- open_expr_tables(count_dir)
#   saveRDS(l_expr_tables, file = "Data/l_expr_tables.rds")
# } else {
#   message ("l_expr_tables has already been created and is sittign in /Data/")
# }
