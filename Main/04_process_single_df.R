#This script processes a single data frame within l_expr_tables. 
# l_expr_tables is opened, and a test_frame is selected from it
# remove gene id versions from test data frame
# merge test data frame with pc table. Name it test_merged
# It then processes the test data frame. This processing will create a df containing expressed genes
# Get some visualization about the filtering is generated.
# Get some summary statistics about our data frame

#---------------------------------------------------------------------------
#Setup and opening files
#---------------------------------------------------------------------------

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

#---Create place to save results
if (!dir.exists("Data/test_frame_analysis")) {
  dir.create("Data/test_frame_analysis")
}

#---Select 1 table from l_expr_tables. Save it as an rds object for later use
#-index is the frame that you are selecting from l_expr_tables
index <- 3
test_frame <- l_expr_tables[[index]]
              
#---Running remove_gene_id_vers on test data frame to remove the gene IDs
#-the expression count tables have gene_id versions, for example, the ".10" in "ENSMUSG00000031965.10"
#-these are artifacts of Ensemble updating the gene versions. They must be removed in order to merch with the protein coding table,
#-which does not have the version ids

test_frame <- remove_gene_id_vers(test_frame)

#---Open pc table

if (file.exists(file = "Data/ensembl_mouse_protein_coding_104.tsv")) {
  message("protein coding table found, opening...")
  pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
} else {
  message("Couldn't find pc table, it should be been downloaded in ")
}

#---------------------------------------------------------------------------
#Perform some minor processing on the protein coding table (pc)
#---------------------------------------------------------------------------

#---Remove nondistinct from pc table and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  select(Gene_ID, Symbol)

#---Identifying gene_ids with more than 1 symbol.
#-The protein coding gene has some 20 Symbols that correspond to more than 1 gene_id. 
#-As the rest of the processing deals focuses on merging via gene_id, and not symbol, these 20 Symbols were left as is associated
#-with their various unique gene ids.
pc_sub_duplicates <- pc_sub %>%
  filter(duplicated(Symbol))%>%
  filter(Symbol !="")%>%
  distinct(Symbol)


#---------------------------------------------------------------------------
#merge pc table onto the expression test frame
#---------------------------------------------------------------------------

test_merged <- test_frame %>%
  dplyr::rename(Gene_ID = gene_id) %>% 
  left_join(y = pc_sub, by = "Gene_ID") %>% 
  relocate(Symbol, .after = Gene_ID)

#---------------------------------------------------------------------------
#Looking at breakdown of protein coding genes that did not join on to the count table
#---------------------------------------------------------------------------
# A quick glance at the protein coding genes that were not expressed in the count table shows that many were
# Olfactory proteins, and RIKEN-named genes
did_not_join <- filter(pc_sub, !pc_sub$Gene_ID %in% test_merged$Gene_ID)

#---------------------------------------------------------------------------
# This part of the  script processes 1 merged expression - pc data frame. Into gene_ids with symbols(pc), and gene_ids without symbols
#---------------------------------------------------------------------------

#---get a data frame of all the empty symbols

empty_symbols <-  filter(test_merged, is.na(Symbol))

#---Filter out transcripts that do have a protein symbol

with_symbols <- filter(test_merged, !(is.na(Symbol)))

#---Check to make sure empty_symbols + with_symbols is equal to the total amount of rows in test_merged
stopifnot(nrow(test_merged) == nrow(empty_symbols)+nrow(with_symbols))

#---Check that all of your transcriopts in with_symbols do indeed have symbols
stopifnot(!(all(is.na(with_symbols$Symbol))))

#---Check that all of your transcriopts in empty_symbols do NOT have symbols
stopifnot((all(is.na(empty_symbols$Symbol))))

#---Make sure the number of distinct genes was equal to the num of groups before
n_groups <- test_merged%>%
  filter(!is.na(Symbol))%>%
  group_by(Gene_ID)%>%
  tally()%>%
  nrow()
stopifnot(nrow(with_symbols) == n_groups)

#---Save with_Symbols as an RDS object
saveRDS(with_symbols, file = "Data/test_frame_analysis/merged_with_Symbols.Rds")



#---------------------------------------------------------------------------
# This part of the script is for generating Alex-requested summary stats
#---------------------------------------------------------------------------

n_matched_to_symbol <- nrow(with_symbols)

frac_matched_to_symbol <- nrow(with_symbols)/nrow(test_merged) 

n_not_matched_to_symbol <- nrow(empty_symbols)

frac_not_matched_to_symbol <-  nrow(empty_symbols)/nrow(test_merged)

stopifnot(frac_matched_to_symbol+frac_not_matched_to_symbol == 1)

n_not_matched_pc_table <- nrow(pc_sub) - sum(pc_sub$Gene_ID %in% test_merged$Gene_ID)

#pc_sub has no empty rows. Each match has a Symbol
#rows removed in n_not_matched_to_symbol did not have a symbol. Likely non-coding genes and such

#---------------------------------------------------------------------------
# This part of the script is for plotting
#---------------------------------------------------------------------------
names <- factor(x = c("Size of PC Table", 
                      "Num Gene Id's that didn't match", 
                      "Num Gene Id's that matched", 
                      "Size of Merged Count Table", 
                      "Num of Expressed Proteins with Symbols", 
                      "Num of Expressed Proteins without Symbols"),
                levels = c("Size of PC Table", 
                           "Num Gene Id's that didn't match", 
                           "Num Gene Id's that matched", 
                           "Size of Merged Count Table", 
                           "Num of Expressed Proteins with Symbols", 
                           "Num of Expressed Proteins without Symbols"))

lengths <- c(nrow(pc_sub), 
             n_not_matched_pc_table, 
             sum(pc_sub$Gene_ID %in% test_merged$Gene_ID),
             nrow(test_merged), 
             n_matched_to_symbol, 
             n_not_matched_to_symbol)

df_summary_stats <- data.frame(x = names, y = lengths)

pl_summary_stats <- ggplot(df_summary_stats) +
  geom_col(mapping = aes(x = names, y = lengths))+
  labs(title = names(l_expr_tables)[[index]])+
  xlab("Summary Stats")+
  ylab("Counts of Gene Id\'s")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
pl_summary_stats
  
write_delim(df_summary_stats, file = paste0("Data/test_frame_analysis/",names(l_expr_tables)[[index]],"df.tsv"), delim = "\t" )
ggsave(file = paste0("Data/test_frame_analysis/",names(l_expr_tables)[[index]],"plot.png"), plot = pl_summary_stats)

#--------------End 












#----------ARCHIVE---------------
# names <- c("Protein Coding Table", "Expression Table", "Merged Frame", "Transcripts with Symbol",
#            "Transcripts with no Symbol", "Unique Expressed Genes")
# names <- factor(names, levels = names)
# 
# lengths <- c(nrow(pc), nrow(test_frame), nrow(test_merged), n_pc_expr_1_Sym, n_pc_expr_1_noSym, nrow(test_merged_uniq))
# 
# summary_frame <- data.frame("names" = names, 
#                             "lengths" = lengths)
# summary_plot <- ggplot(summary_frame, mapping = aes(x = names, y = lengths)) +
#   geom_col()
# 
# #---Save plot as rds object
# 
# if (!dir.exists("Plots")) {
#   dir.create("Plots")
# }
# 
# saveRDS(summary_plot, file = "Plots/summary_plot")


# #---Get summary statistics for test_merged_uniq
# 
# summary_stats <- test_merged_uniq %>%
#   select(Symbol, expected_count, TPM, FPKM)%>%
#   skim()
# summary_stats
# 
# expected_count_dist <- ggplot(test_merged_uniq)+
#   geom_histogram(mapping = aes(expected_count), binwidth = 1000)+
#   scale_y_continuous(trans='log10')
# TPM_dist <- ggplot(test_merged_uniq)+
#   geom_histogram(mapping = aes(TPM))+
#   scale_y_log10()
# FPKM_dist <- ggplot(test_merged_uniq)+
#   geom_histogram(mapping = aes(FPKM))+
#   scale_y_log10()
# test_merged_uniq$FPKM
# 
# ggplot(data = test_merged_uniq, mapping = aes(expected_count)) +
#   geom_boxplot()+
#   scale_x_continuous(trans = 'log10')
# ggplot(data = test_merged_uniq, mapping = aes(TPM)) +
#   geom_boxplot()+
#   scale_x_continuous(trans = 'log10')
# ggplot(data = test_merged_uniq, mapping = aes(FPKM)) +
#   geom_boxplot()+
#   scale_x_continuous(trans = 'log10')





# #add this later on to the iterative version of remove gen id version
# 
# #
# #Check to ensure that remove_gene_id_vers was successful for all tables
# #
# stopifnot(!FALSE %in% (lapply(l_expr_tables,check_gene_id_vers_removal)))
