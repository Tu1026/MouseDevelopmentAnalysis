#This script processes a single data frame within l_expr_tables. 
# l_expr_tables is opened, and a test_frame is selected from it
# remove gene id versions from test data frame
# merge test data frame with pc table. Name it test_merged
# It then processes the test data frame. This processing will create a df containing expressed genes
# Get some visualization about the filtering is generated.
# Get some summary statistics about our data frame

# Returns:
#test_merged. The expression table test merged onto pc table
#test_merged_unique. The expression table containing only unique expressed genes
#summary_visualization. visualization plot 


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
test_frame <- l_expr_tables[[1]]
              
if(!file.exists("Data/test_frame_analysis/test_frame.tsv")) {
  file = write_tsv(test_frame, file = "Data/test_frame_analysis/test_frame.tsv")
}


#---Running remove_gene_id_vers on test data frame to remove the gene IDs

test_frame <- remove_gene_id_vers(test_frame)

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
  # filter(Symbol == "Pcdha11" )
# nrow(pc_sub)
# # glimpse(pc_sub)
# # ok <- filter(pc_sub, pc_sub$Symbol =="Pcdha11" )
# # distinct(ok, Symbol)
# # pc_sub_uniq <- unique(pc$Symbol)
# # length(pc_sub_uniq)
#----These 2 are off by 48 rows
# nrow(pc_sub)-length(pc_sub_uniq)

# Here we can see the 48 rows. We can see that Pakap and pcdha11 are our culprits here
# glimpse(pc_sub[duplicated(pc_sub$Symbol),])


#---merge pc table with test data frame

test_merged <- test_frame %>%
  dplyr::rename(Gene_ID = gene_id) %>% 
  left_join(y = pc_sub, by = "Gene_ID") %>% 
  relocate(Symbol, .after = Gene_ID)

#---save test_merged as a tsv
if(!file.exists("Data/test_frame_analysis/test_merged.tsv")) {
  write_tsv(x = test_merged, file = "Data/test_frame_analysis/test_merged.tsv")
}


#---------------------------------------------------------------------------
# This part of the  script processes 1 merged expression - pc data frame. 
#---------------------------------------------------------------------------

#---get a data frame of all the empty symbols

empty_symbols <- test_merged%>%
  filter(is.na(Symbol))

#---Filter out transcripts that do have a protein symbol

with_symbols <- test_merged %>%
  filter(!(is.na(Symbol)))

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

# test_frame[duplicated(test_frame$gene_id),]
# test_merged[duplicated(test_merged$Symbol),]
# filter(test_merged, test_merged$Symbol=="Pcdha11")
# test_with_symbols <- test_merged %>%
#   filter(!(is.na(Symbol)))
# glimpse(test_with_symbols[duplicated(test_with_symbols$Symbol),])
# glimpse(pc_sub[duplicated(pc_sub$Symbol),])

#---Save with_Symbols as an RDS object
saveRDS(with_symbols, file = "Data/test_frame_analysis/merged_with_Symbols.Rds")

#---------------------------------------------------------------------------
# This part of the script is for generating Alex-requested summary stats
#---------------------------------------------------------------------------

n_matched_to_symbol <- with_symbols %>%
  nrow()
frac_matched_to_symbol <- with_symbols %>%
  nrow()/nrow(test_merged)
n_not_matched_to_symbol <- empty_symbols %>%
  nrow()
frac_not_matched_to_symbol <- empty_symbols %>%
  nrow()/nrow(test_merged)
stopifnot(frac_matched_to_symbol+frac_not_matched_to_symbol == 1)

n_not_matched_pc_table <- nrow(pc_sub) - sum(pc_sub$Gene_ID %in% test_merged$Gene_ID)

#pc_sub has no empty rows. Each match has a Symbol
#rows removed in n_not_matched_to_symbol did not have a symbol. Likely non-coding genes and such

#---------------------------------------------------------------------------
# This part of the script is for plotting
#---------------------------------------------------------------------------
names <- factor(x = c("PC", "Gene Id's didn't match", "Gene Id's matched", "Merged Count Table", "Merged with Symbols", "Merged without Symbols"),
                levels = c("PC", "Gene Id's didn't match", "Gene Id's matched", "Merged Count Table", "Merged with Symbols", "Merged without Symbols"))
lengths <- c(nrow(pc_sub), n_not_matched_pc_table, sum(pc_sub$Gene_ID %in% test_merged$Gene_ID), nrow(test_merged), n_matched_to_symbol, 
             n_not_matched_to_symbol)
plot_frame <- data.frame(x = names, y = lengths)
ggplot(plot_frame) +
  geom_col(mapping = aes(x = names, y = lengths))


#--------------End 












#----------ARCHIVE---------------
names <- c("Protein Coding Table", "Expression Table", "Merged Frame", "Transcripts with Symbol",
           "Transcripts with no Symbol", "Unique Expressed Genes")
names <- factor(names, levels = names)

lengths <- c(nrow(pc), nrow(test_frame), nrow(test_merged), n_pc_expr_1_Sym, n_pc_expr_1_noSym, nrow(test_merged_uniq))

summary_frame <- data.frame("names" = names, 
                            "lengths" = lengths)
summary_plot <- ggplot(summary_frame, mapping = aes(x = names, y = lengths)) +
  geom_col()

#---Save plot as rds object

if (!dir.exists("Plots")) {
  dir.create("Plots")
}

saveRDS(summary_plot, file = "Plots/summary_plot")


#---Get summary statistics for test_merged_uniq

summary_stats <- test_merged_uniq %>%
  select(Symbol, expected_count, TPM, FPKM)%>%
  skim()
summary_stats

expected_count_dist <- ggplot(test_merged_uniq)+
  geom_histogram(mapping = aes(expected_count), binwidth = 1000)+
  scale_y_continuous(trans='log10')
TPM_dist <- ggplot(test_merged_uniq)+
  geom_histogram(mapping = aes(TPM))+
  scale_y_log10()
FPKM_dist <- ggplot(test_merged_uniq)+
  geom_histogram(mapping = aes(FPKM))+
  scale_y_log10()
test_merged_uniq$FPKM

ggplot(data = test_merged_uniq, mapping = aes(expected_count)) +
  geom_boxplot()+
  scale_x_continuous(trans = 'log10')
ggplot(data = test_merged_uniq, mapping = aes(TPM)) +
  geom_boxplot()+
  scale_x_continuous(trans = 'log10')
ggplot(data = test_merged_uniq, mapping = aes(FPKM)) +
  geom_boxplot()+
  scale_x_continuous(trans = 'log10')

#TODO Need advice on how to visualize




# #add this later on to the iterative version of remove gen id version
# 
# #
# #Check to ensure that remove_gene_id_vers was successful for all tables
# #
# stopifnot(!FALSE %in% (lapply(l_expr_tables,check_gene_id_vers_removal)))
