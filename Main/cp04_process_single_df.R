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
tryCatch(
  {
    l_expr_tables <- readRDS(file = "Data/l_expr_tables")
  },
  error = function(cond) {
    message(cond)
    message("Couldn't find l_expr_tables rds object. Should have been generated in 01")
  }
)

#---Select 1 table from l_expr_tables. Save it as an rds object for later use

test_frame <- l_expr_tables[[1]]

saveRDS(test_frame, file = "Data/test_frame")
#---Running remove_gene_id_vers on test data frame to remove the gene IDs

test_frame <- remove_gene_id_vers(test_frame)

#---Open pc table
tryCatch(
  {
    pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
  },
  error = function(cond) {
    message(cond)
    message("Couldn't open the pc table. It should be in Data/")
  }
)


#---merge pc table with test data frame

# test_merged <- left_join(pc, test_frame, by= c("Gene_ID" = "gene_id"))

pc_sub <- pc %>% 
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  select(Gene_ID, Symbol)
  
test_merged <- test_frame %>%
  dplyr::rename(Gene_ID = gene_id) %>% 
  left_join(y = pc_sub, by = "Gene_ID") %>% 
  relocate(Symbol, .after = Gene_ID)


#---save test_merged as rds object

saveRDS(test_merged, file = "Data/test_merged")

#-----------------------------
# This part of the  script processes 1 merged expression - pc data frame. 

#---Check to see if there are any empty symbols

stopifnot(all(!(is.na(test_merged$Symbol))))

#---Filter out transcripts that do have a protein symbol

pc_expr_1_Sym <- test_merged %>%
  filter(!(is.na(Symbol)))

#---Get the number of transcripts with symbols
n_pc_expr_1_Sym <- pc_expr_1_Sym %>%
  nrow()

#---Check that all of your transcriopts do indeed have symbols
stopifnot(!(all(is.na(pc_expr_1_Sym$Symbol))))


#---Filter out transcripts that NOT have a protein symbol
pc_expr_1_noSym <- test_merged %>%
  filter((is.na(Symbol)))

#---Get the number of transcripts with symbols
n_pc_expr_1_noSym <- pc_expr_1_noSym %>%
  nrow()

#---Check that all of your transcriopts do NOT have symbols
stopifnot((all(is_empty(pc_expr_1_noSym$Symbol))))

#---Check to make sure your number of symbol transcripts + non-symbol transcripts adds up to total
stopifnot(n_pc_expr_1_noSym + n_pc_expr_1_Sym == nrow(test_merged))


#---get only the unique symbols. Delete the duplicates
test_merged_uniq <- test_merged %>%
  distinct(Symbol, .keep_all = TRUE)

#---Make sure the number of distinct genes was equal to the num of groups before
n_groups <- test_merged%>%
  group_by(Symbol)%>%
  tally()%>%
  nrow()
stopifnot(nrow(test_merged_uniq) == n_groups)

#---Save test_merged_uniq as an RDS object
saveRDS(test_merged_uniq, file = "Data/unique_test_merged")

#---Basic plotting to get an idea of how the poportion of matches we have to the total pc table

#--first we have to open up pc table again

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
