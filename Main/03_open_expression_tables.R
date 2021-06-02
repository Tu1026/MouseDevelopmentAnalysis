#
#
source('Main/libraries.R')

#This script opens up l_expr_tables. A list containing all of the downloaded expression tables.
#This script also removes the gene_id transcript versions from all expression tables

if (!"l_expr_tables" %in% ls()) {
 return(l_expr_tables <- open_expr_tables(count_dir))
} else {
  warning ("l_expr_tables has already been created")
}


#
#Removal Events
#
