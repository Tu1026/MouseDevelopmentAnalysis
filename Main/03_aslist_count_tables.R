#This script opens all of the expression tables in Data/count_tables as a list.

#Returns:
#It then saves this list as an rds object named l_expr_tables. 


source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(stringr)


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
# ordered in the same order as the meta_data table. This is for convenience


expr_table_names <- list() 
for (i in meta_data$id) {
  i <- paste0(i,".tsv")
  expr_table_names <- append(expr_table_names, i)
}

#---------------------------------------------------------------------------
# Open the count tables
#--------------------------------------------------------------------------
# Open count tables and save inside a list. List should be ordered
# In the order of the meta_data table.
# Remove the .tsv after the names once the files are read

l_expr_tables <- lapply(expr_table_names, open_expr_table, target_directory)

expr_table_names <- meta_data$id
names(l_expr_tables) <- expr_table_names

saveRDS(l_expr_tables, file = "Data/l_expr_tables.rds")


# #---Rename list intexes to table names
# rename_expr_tables <- function(l_expr_tables, list_of_expression_tables) {
#   names <- str_replace(list_of_expression_tables, pattern = ".tsv", replacement = "")
#   names(l_expr_tables) <- names
#   return(l_expr_tables)
# }

##---function to open all tables within 1 list
# open_expr_tables <- function(target_directory, list_of_expression_tables) {
#   #Open up all of the expression tables in list_of_expression_tables.
#   #all tables must be tsvs, and they must all be in target_directory
#   #return l_expr_tables, which is also saved as an rds in target_directory
#   
#   l_expr_tables <- lapply(list_of_expression_tables, open_count_table, target_directory = target_directory)
#   l_expr_tables <- rename_count_tables(l_expr_tables, list_of_expression_tables)
#   return(l_expr_tables)
# }





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
