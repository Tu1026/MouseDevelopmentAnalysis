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
#rename pc table's Gene_ID to gene_id for future use
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
# Create l_pc_expr_tables. A list containing tables of the expression of pc genes
#

#We want to join the pc table to our expression table. This will result in duplicate lines for each of the duplicate gene_ids in the pc table.
#But it will avoid the creation of a bunch of nulls for values in the pc table

create_pc_expr_tables <- function(pc_table,expression_df) {
  #
  #joins the protein coding table onto the expression coding table by gene_id.
  #
  #@param: pc_table: the protein coding table downlaoded by ensemble. Should be named "pc"
  #@param: expression_df: A given expression data frame. Should be within "l_expr_tables"
  #return: new data frame thats merged both tables merged by "gene_id" columns
  stopifnot("gene_id" %in% names(pc_table) & "gene_id" %in% names(expression_df))
  merged_tables <- right_join(pc_table,expression_df, by = "gene_id")
  return (merged_tables)
  }

l_pc_expr <- lapply(l_expr_tables, create_pc_expr_tables, pc_table=pc)


glimpse(l_expr_tables[[1]])
glimpse(l_pc_expr[[1]])



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
# 
# 
# 
# #
# #Compare expression tables to ensemble's PC. Return only PC genes
# #
# 
# glimpse(pc)
# 
# 
# 
# # expr_dat <- read.delim(file = paste0(count_dir,"/","ENCFF227HKF.tsv"))
# # gene_ids <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
# # num_pc <- table(gene_ids %in% pc$Gene_ID)
# # pc_gene_ex <- expr_dat %>%
# #   select(c(gene_id,length,expected_count))%>%
# #   mutate(gene_id = gene_ids) %>%
# #   rename(Gene_ID = gene_id) %>%
# #   filter(Gene_ID %in% pc$Gene_ID)%>%
# #   left_join(y=pc, by="Gene_ID")%>%
# #   select(c(Gene_ID,length,expected_count,Chromosome,Symbol))%>%
# #   distinct()
# # glimpse(pc_gene_ex)  
# 
# 
# # Overloaded function with no documentation. Consider function that acts on
# # 'unit.' Then some other iterative code calls that function for every unit
# # Also way too prelim. Need to take time to inspect presence of genes for each
# # table. Do they all have equal genes? How many are being thrown away? What 
# # if one PC is present in one table but not another? Why did you grab expected
# # counts with no mention of TPM/RPKM? You call distinct, which would just blindly
# # throw away some transcripts (would keep the first element in encounters - is
# # this desirable?). Avoid functions that are writing data without being
# # explicit about it
# 
# 
# get_pc_genes <- function(gen_exp_files_names) {
#   if(!(dir.exists("Data/Processed_Count_tables"))){
#     warning("Creating Data/Processed_Count_tables directory")
#     dir.create("Data/Processed_Count_tables")
#   }
#   for (i in gen_exp_files_names) { 
#     expr_dat <- read.delim(file = paste0(count_dir,"/",i))
#     gene_ids <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
#     num_pc <- table(gene_ids %in% pc$Gene_ID)
#     pc_gene_ex <- expr_dat %>%
#       select(c(gene_id,length,expected_count))%>%
#       mutate(gene_id = gene_ids) %>%
#       rename(Gene_ID = gene_id) %>%
#       filter(Gene_ID %in% pc$Gene_ID)%>%
#       left_join(y=pc, by="Gene_ID")%>%
#       select(c(Gene_ID,length,expected_count,Chromosome,Symbol))%>%
#       distinct() %>%
#       mutate(name = i)
#     # stopifnot(gene_ids == pc_gene_ex$Gene_ID)
#     # stopifnot(nrow(pc_gene_ex)==num_pc[2])
#     # glimpse(pc_gene_ex)
#     write_tsv(x=pc_gene_ex,
#               file=paste0("Data","/","Processed_Count_tables","/","pc_gene_ex_",as.character(i)))
#   }
# }
# my_expression_tables <- list.files(count_dir)
# 
# # AM: This check doesn't do anything. A string is not null
# 
# if (is.null("Data/Processed_Count_tables")){
#   get_pc_genes(my_expression_tables)
# }
# 
# #----------  open all processed count data
# 
# # open_pro_tables <- function() {
# #   tables <- c()
# #   for (i in my_expression_tables) {
# #     print (i)
# #     name <- paste0("pro_",i)
# #     table <- read.delim(paste0("Data/Processed_Count_tables/pc_gene_ex_",i))
# #     print(name)
# #     glimpse(table)
# #     name <-  table
# #     append(x=tables, values=name)
# #   return(tables)
# #   }
# # }
# # 
# # open_pro_tables()
# 
# pro_dir <- "Data/Processed_Count_tables"
# 
# processed_count_tables <- paste0(pro_dir,"/",list.files(pro_dir))
# list_pro_count_tables <- lapply(processed_count_tables,read.delim)
# names(list_pro_count_tables) <- list.files(pro_dir)
# 
# #todo make this more robust with regex
# stopifnot (length(list_pro_count_tables)== length(my_expression_tables))
# 
# 
# #---------- apply replicate metadata from google sheets
# 
# # AM: Doesn't work - need to load appropriate library
# # Also naming needs to include E/P
# # Don't need to add metadata into expression table: keep separate metadata
# 
# replicate_meta <- read_sheet("1gFkSHD15wdd3FdCrZ51DQOi9RYQg01i1TXpIvkDVIz0")
# 
# 
# correct_levels <- c(10.5,11.5,12.5,13.5,14.5,15.5,16.5,0)
# 
# identify_in <- function(dataframe,replicate_meta){
#   my_row <- replicate_meta[which(replicate_meta$file_dataset == dataframe$name[1]),]
#   dataframe <- dataframe%>%
#     mutate(tissue_type = factor(my_row[2]),
#            dev_stage = factor(my_row[3], levels = correct_levels),
#            replicate = factor(my_row[4]))
# }
# 
# list_pro_count_tables <- lapply(list_pro_count_tables,identify_in,replicate_meta)
# # list_pro_count_tables
# # glimpse(list_pro_count_tables[1])
# # head(list_pro_count_tables[1])
# 
# 
# # replicate_meta[which(replicate_meta$file_dataset == "ENCFF227HKF.tsv"),]
# # list_pro_count_tables[names(list_pro_count_tables)[1]]
# 
# 
# 
# #----------- conglomerate important info into 1 master table
# # create_master_table <- function(list_of_dataframes) {
# #   master_table <- data.frame()
# #   for (i in list_of_dataframes){
# #     print(i)
# #     # print(master_table)
# #     # master_table <- rbind(master_table,i)
# #   return(master_table)
# #   }
# # }
# 
# #why was that not working????
# 
# master_table <- data.frame()
# for (i in list_pro_count_tables){
#   # print(i)
#   # print(master_table)
#   master_table <- rbind(master_table,i)
#   # print(master_table)
# }
# 
# 
# nrows <-lapply(list_pro_count_tables, function(list_of_tables){nrow(list_of_tables)})
# count_table_rows <- sum(unlist(nrows))
# 
# stopifnot(nrow(master_table)==count_table_rows)
# 
# #---------- summary stats on master_table
# 
# #todo pearson 
# 
# sum_master_table <- master_table %>%
#   group_by(dev_stage,Symbol)%>%
#   mutate(replicate=as.double(replicate))%>%
#   summarize(avg_count = mean(expected_count), pearson = cor(replicate,expected_count,method="pearson"))
# 
# glimpse(sum_master_table)
# 
# arr_master_table <- sum_master_table %>%
#   arrange(desc(avg_count), .by_group=TRUE)
# 
# 
# #------------prelim visualization
# 
# library(ggplot2)
# 
# sample_table <- arr_master_table %>%
#   filter(dev_stage == c(10.5,11.5),
#          avg_count > 10000)
# glimpse(sample_table)
# 
# dev_stage_graph <- ggplot(sample_table)+
#   geom_col(mapping = aes(x = Symbol, y = avg_count, fill = dev_stage), position = 'dodge')+
#   theme(axis.ticks.x= element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank())