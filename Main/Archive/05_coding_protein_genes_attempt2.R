#
#01, and 03 must be run prior to this script
#

#This script takes l_expr_tables from 03, and uses pc from 01 to filter for protein coding genes


library(tidyverse)

glimpse(pc)
glimpse(l_expr_tables)

test_dataframe <- l_expr_tables[[1]]
glimpse(test_dataframe)
glimpse(test_dataframe %>%
          dplyr::filter((test_dataframe$gene_id == test_dataframe$transcript_id.s.) == FALSE))
#
#rename pc table's Gene_ID to gene_id for future use for merging
#
if ("Gene_ID" %in% names(pc)) {
  pc <- rename(pc,"gene_id"="Gene_ID")
}
stopifnot(names(pc[6]) == "gene_id")

#
#Identifying how many of our transcripts in our expression are pc genes
#
summarize_genes <- function(expr_df, pc_table) {
  
  #
  #Creates a sumamry table of information about an expresion table
  #
  #@param: expr_dfL An expression data frame. Must have its transcript versions removed and have a gene_id column
  #@param: pc_table: protein coding table. Must have its Gene_ID col changed to gene_id
  
  tot_transcripts <-  nrow(expr_df)
  all_pc_genes <- filter(expr_df, str_detect(expr_df$gene_id, "ENSM"))
  tot_pc_genes <-  nrow(filter(expr_df, str_detect(expr_df$gene_id, "ENSM")))
  tot_pc_genes_in_expr <-  nrow(filter(expr_df, expr_df$gene_id %in% pc_table$gene_id))
  tot_pc_genes_not_in_expr <- nrow(filter(all_pc_genes, !(all_pc_genes$gene_id %in% pc_table$gene_id)))
  tot_not_pc_gene_transcripts <- tot_transcripts - tot_pc_genes
  return (data.frame(tot_transcripts, tot_pc_genes, tot_pc_genes_in_expr, tot_pc_genes_not_in_expr, tot_not_pc_gene_transcripts))
}

l_expr_sum <- lapply(l_expr_tables, summarize_genes , pc_table = pc)

l_expr_sum
l_expr_sum[[1]]

#
#Check to make sure that the tot_genes is equal to tot_pc genes and tot_not_pc_genes
#

check_tot_gene <- function(expr_sum_df) {
  #
  #checks to make sure that the tot_genes == tot_pc+tot_not_pc_genes
  #
  
  #@param: A single data frame containing summarized information about an expression data frame. Basically a df in l_expr_sum
  return (stopifnot(expr_sum_df$tot_pc_genes == (expr_sum_df$tot_pc_genes_in_expr + expr_sum_df$tot_pc_genes_not_in_expr)))
}

check_tot <- lapply(l_expr_sum, check_tot_gene)
stopifnot(all(for (i in check_tot) {
  is.null(i)
}))



#It looks like each expression table has 21302 protein coding genes.



#
# l_pc_expr has data frames with the pc merged onto the expression data frames. But it seems like there are gene_ids on the expression data frames
# that are actually in the expression data frames. But when the pc is merged on, these weird gene_id rows are not deleted. They just have NA on them
#the entire way. This is just to see how many of these gene_ids are present. We will do this by filtering by symbol==na

my_df <- l_pc_expr[[1]]

actual_pc_genes <- my_df %>%
  filter(!is.na(my_df$Symbol))
non_pc_genes <- my_df %>%
  filter(is.na(my_df$Symbol))

num_actual_pc_genes <- actual_pc_genes %>%
  nrow()
num_non_pc_genes <- non_pc_genes %>%
  nrow()

#
#Now that I have the number, I'm going to filter out all of the non_pc_genes
#

filter_non_pc <- function(data_frame) {
  #
  # Purpose: Filter out rows that have Symbol == NA. These all are not real protein coding genes
  #
  #@param: data_frame that is in l_pc_expr
  #@return: data frame with rows removed where there is no protein 
  data_frame <- data_frame %>%
    filter (!is.na(data_frame$Symbol))
  return (data_frame)
}

#
# Executing that filter_non_pc to all data frames in the protein coding data frame
#

l_pc_expr_fil <- lapply(l_pc_expr,filter_non_pc)

#
# Now lets check
#

my_test <- l_pc_expr_fil[[1]]
glimpse(my_test)

my_test %>%
  group_by(Symbol)

tot_pc_genes_in_expr <-  filter(test_dataframe,test_dataframe$gene_id %in% pc$gene_id)

please_work <- filter(my_test, (my_test$gene_id %in% tot_pc_genes_in_expr$gene_id))
help_me <- filter(tot_pc_genes_in_expr, tot_pc_genes_in_expr$gene_id %in% my_test$gene_id)
nrow(help_me)
please_work <- group_by(please_work, Symbol)
nrow(please_work)

nrow(tot_pc_genes_in_expr)


#
#------ each of our protein coding tables contains only protein coding genes theoretically.
#
#But it seems like there are some weird gene_ids in the pc table. Such as "9641". It doesn't seem to be a pc gene.
#These values are present in our expression tables, but their counts are always 0. Take a look at below
pc_expr <- l_pc_expr[[1]]
not_ens <- (pc_expr %>%
              filter(!str_detect(pc_expr$gene_id, "ENS") == TRUE))
is_ens <- (pc_expr %>%
             filter(str_detect(pc_expr$gene_id, "ENS") == TRUE))
is_ens$gene_id


nrow(pc_expr)
nrow(not_ens)
nrow(is_ens)
stopifnot(nrow(pc_expr)==nrow(not_ens)+nrow(is_ens))

#from not_ens, you can see that there are a bunch of weird "pc genes" that dont have ENS values as gene_ids
not_ens$gene_id
#you can see that these values are indeed in the pc table
all(not_ens$gene_id %in% pc_expr$gene_id)
#You can also see that these gene's have 0 expression count in the expression table
proof <- pc_expr%>%
  filter(pc_expr$gene_id %in% not_ens$gene_id)
proof$expected_count

all(not_ens$gene_id %in% l_pc_expr[[1]]$gene_id)


#the fact that this turns TRUE, it means that the there's all of these l_expr_tables were merged with pc. But for some reason, the l_pc_expr frames shouldn't have
#any genes with no symbol. But there are some, such as gene id "11094"
all(not_ens$gene_id %in% l_pc_expr[[1]]$gene_id)



####################* Because in summarize_genes(), I use str_detect with pattern "ENS" to detect protein coding genes. 
####################*These genes would not be accounted for as protein coding genes. This is correct I believe
####################*However, it must be noted that before we confirm that our expression tables have 21302 pc genes, we need to first delete
####################*All of those fake genes that don't start with ENS

#
#Deleting non-ENS gene_ids from expression tables
#

delete_non_ens <- function(df_pc_expr) { 
  #
  #delete all gene_ids that don't start with "ENS"
  #
  not_ens <- (df_pc_expr %>%
                filter(!str_detect(df_pc_expr$gene_id, "ENS") == TRUE))
  is_ens <- (df_pc_expr %>%
               filter(str_detect(df_pc_expr$gene_id, "ENS") == TRUE))
  stopifnot(nrow(df_pc_expr) == nrow(not_ens)+ nrow(is_ens))
  return (is_ens)
}

l_pc_expr <- lapply(l_pc_expr, delete_non_ens)

glimpse(l_pc_expr[[1]])
#now, all of the expression tables should be linked to pc info. And the total amount of unique proteins should be 21302. 
#If not I'm going to be sad

#
# Checking count of proten coding genes in expression data
#
my_pc_expr1 <- l_pc_expr[[1]]%>%
  group_by(Symbol)
n_uniq_pc <- n(pc_table)
glimpse(pc_table)

filter(pc_table, expr_df$gene_id %in% pc_table$gene_id)
all(pc_table$gene_id )


#not sure what these pc genes are. I guess I'll have to consult ensemble. Regardless, they should likely be deleted so I'll delete them from the expr tables

#
# Deleting the weird non-pc genes
#


# #
# #filtering out non-pc genes
# #
# 
# 
# 
# l_pc_expr[[1]]%>%
#   group_by(Symbol)
# 
# n <- 21786
# 
# 
# l_pc_in_expr[1]
# #TODO Why is this not working? my n() is returning fatal error
# 
# #
# # Taking a look at the merged tables
# #
# 
# 