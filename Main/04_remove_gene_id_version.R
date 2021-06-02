#
#
library(tidyverse)

###This script removes the decimals after each value of the expression table's gene_id column

# # TODO: What exactly do these diff versions mean and should we remove them?


#
#Running remove_gene_id_vers on all expression data dables
#

l_expr_tables <- lapply(l_expr_tables,remove_gene_id_vers)

#
#Check to ensure that remove_gene_id_vers was successful for all tables
#
stopifnot(!FALSE %in% (lapply(l_expr_tables,check_gene_id_vers_removal)))



#
#Removal Events
#
