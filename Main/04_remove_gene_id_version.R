#
#requires 01, 03 to be run
#

whats_run <- c(whats_run,"04")
library(tidyverse)

###This script removes the decimals after each value of the expression table's gene_id column


#######IGNORE THIS COMMENTED BLOCK

# ###It looks like each ENS has a geneId transcript version. there are no ENSMUSG[:digits:], all are ENSMUSG[:digits:].3 or .2 or .6 etc
# l1 <- length(str_subset(test_ENCFF227HKF.tsv$gene_id,"ENS"))
# l2 <- length(str_subset(test_ENCFF227HKF.tsv$gene_id,"\\."))
# stopifnot(l1 == l2)
#
#
# #

#
# ###############################CHECK TO SEE IF FUNCTION WORKS WITH A N=50 size
# ###Get Sample to check with
# my_pre_test_sample <- sample_n(test_ENCFF227HKF.tsv, size = 50)
# my_pre_test_sample$gene_id
# str_subset(my_pre_test_sample$gene_id, "ENS")
# my_pre_test_values <- c("ENSMUSG00000014418.11","ENSMUSG00000045395.3","ENSMUSG00000084007.3","ENSMUSG00000016758.3")
# my_post_test_values <- c("ENSMUSG00000014418","ENSMUSG00000045395","ENSMUSG00000084007","ENSMUSG00000016758")
# 
# ###Positive Check Pre-Test: Should be TRUE
# all(my_pre_test_values %in% my_pre_test_sample$gene_id)
# ###Negative Check Pre-Test: Should be false
# all(my_post_test_values %in% my_pre_test_sample$gene_id)
# 
# ###Execution of function
# my_post_test_sample <- remove_gene_id_vers(my_pre_test_sample)
# my_post_test_sample$gene_id
# 
# ###Check to see if the positive check is still true: Should now be false
# all(my_pre_test_values %in% my_post_test_sample$gene_id)
# ###Check to see if the negative check is still false: Should now be False
# all(my_post_test_values %in% my_post_test_sample$gene_id)



#------------------remove different versions of gene_id in the tables within l_expr_tables
# ##TODO: What exactly do these diff versions mean and should we remove them?
#
#
remove_gene_id_vers <- function(data_frame) {
  #
  #Removes the version ids of the gene_id column within all expr tables
  #
  #@param : data frame for expression
  #@return : data frame with gene_id versions removes
  data_frame$gene_id <- str_replace(data_frame$gene_id, "\\.[:digit:]*", "")
  return(data_frame)
}

#
#Generating Positive and Negative Checks
#
check_gene_id_vers_removal <- function(dataframe) {
  #
  #Checks to see that remove_gene_id_vers was successful
  #
  #@param: Dataframe
  #@return: TRUE= Removal was successfull. FALSE= removal failed
  my_checks <- c( "ENSMUSG00000021028.7  ENSMUSG00000032667.13 ENSMUSG00000016758.3  ENSMUSG00000014418.11 ENSMUSG00000092505.1  ENSMUSG00000045395.3 
                ENSMUSG00000084007.3  ENSMUSG00000087052.1  ENSMUSG00000096615.1  ENSMUSG00000081432.1  ENSMUSG00000076888.1  ENSMUSG00000071816.5 
                ENSMUSG00000064860.1  ENSMUSG00000062098.7  ENSMUSG00000093867.3  ENSMUSG00000027257.9  ENSMUSG00000093806.1  ENSMUSG00000087496.1 
                ENSMUSG00000089999.1  ENSMUSG00000029072.3  ENSMUSG00000101458.1  ENSMUSG00000025747.8  ENSMUSG00000095015.3  ENSMUSG00000090383.1 
                ENSMUSG00000057461.1  ENSMUSG00000081266.1  ENSMUSG00000076485.2")
  my_negative_checks <- str_extract_all(my_checks, boundary("word")) #After running remove_gene_id_vers on a data table, this check should yield false
  my_positive_checks <- str_replace(my_negative_checks[[1]], "\\.[:digit:]*", "") #After running remove_gene_id_vers on a data table, this check should yield TRUE
  if (all(!my_negative_checks %in% dataframe$gene_id) == TRUE & all(my_positive_checks %in% dataframe$gene_id) == TRUE) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}


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
rm("remove_gene_id_vers")
