#
#Requires that 01_downloaded_gene_coount_tables.R is run
#

#This script opens up l_expr_tables. A list containing all of the downloaded expression tables.
#This script also removes the gene_id transcript versions from all expression tables


library(tidyverse)
library(tools)

whats_run <- c(whats_run,"03")


#
#lapply to create a list of all expression tables. l_expr_tables is a list containing all expression tables
#

create_names_of_tables<- function(count_dir) {
  #
  #Create a list of all the names of the tables within count_dir
  #
  #@param: directory where expression tables are
  #@return: list of all the names of your expression tables. Includes location in directory
  my_list <- list()
  for (i in list.files(count_dir)) { 
    my_list <- append(my_list, paste0(count_dir,"/",i))
  }
  return (my_list)
}

#
#Open a list of tables using the return of create_list_of_tables
#

open_expr_tables <- function(count_dir) {
  #
  #creates a list containg all expression table dataframes
  #
  #@param: String of the directory name where data frames live. This is passed into create_names_of_tables.
  #example: count_dir <- "Data/Count_tables"
  #@return: A list containing all of the expression table dataframes. They are named according to their name within their directory
  names_of_tables <- create_names_of_tables(count_dir)
  warning("Opening Tables...")
  l_expr_tables <- lapply(names_of_tables,read.delim)
  names(l_expr_tables) <-  list.files(count_dir)
  return (l_expr_tables)
}

if (!"l_expr_tables" %in% ls()) {
 return(l_expr_tables <- open_expr_tables(count_dir))
} else {
  warning ("l_expr_tables has already been created")
}


#
#Removal Events
#
rm("open_expr_tables","create_names_of_tables")
