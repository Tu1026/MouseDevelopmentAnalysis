#This script does the following
#1) Open all of the count tables in Data/count_tables as a list named l_expr_tables
#2) Processes a single count table
#3) Generate Count-sample matrix

source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(stringr)
library(skimr)
library(Hmisc)
library(corrplot)


#---------------------------------------------------------------------------
########################## 1) Opening Count Tables
#--------------------------------------------------------------------------

# Opens all of the count tables in Data/count_tables as a list named l_expr_tables

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
# ordered in the same order as the meta_data table. This is to 
# represent the samples in order of develempental times

ordered_meta_data_id_names <- paste0(meta_data$id,".tsv")

#---------------------------------------------------------------------------
# Open the count tables
#--------------------------------------------------------------------------
# Open count tables and save inside a list. List should be ordered
# In the order of the meta_data table.
# Remove the .tsv after the names once the files are read

l_expr_tables <- lapply(ordered_meta_data_id_names,
                        open_expr_table,
                        target_directory)

names(l_expr_tables) <- str_replace(ordered_meta_data_id_names, ".tsv","")



#---------------------------------------------------------------------------
########################## 2) Processing 1 count table
#--------------------------------------------------------------------------

# Process a single data frame within l_expr_tables. 
# Remove gene id versions from test data frame
# Merge test data frame with pc table. Name it test_merged
# It then processes the test data frame. This processing will create a 
# df containing Measured genes
# Get some visualization about the filtering is generated.
# Get some summary statistics about our data frame

#---------------------------------------------------------------------------
#Setup and opening files
#---------------------------------------------------------------------------

#---Choose which count table to process
# index is the frame that you are selecting from l_expr_tables
index <- 3
test_table <- l_expr_tables[[index]]

#---Create place to save results
if (!dir.exists("Data/test_table_analysis")) {
  dir.create("Data/test_table_analysis")
}

#---Open pc table
pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)

#---------------------------------------------------------------------------
# Running remove_gene_id_vers on test data frame to remove the gene IDs
#---------------------------------------------------------------------------
#-the expression count tables have gene_id versions, 
# for example, the ".10" in "ENSMUSG00000031965.10"
#-These exist due to Ensemble updating the gene versions.
#-They must be removed in order to merge with the protein coding table
# via Gene_ID, which does not have the version ids

test_table <- remove_gene_id_vers(test_table)

#---------------------------------------------------------------------------
#Perform some minor processing on the protein coding table (pc)
#---------------------------------------------------------------------------

#---Remove nondistinct Symbols and Gene_IDs from pc table
# and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)

#---Why we also remove duplicate Symbols...
#-The protein coding gene has some 20 Symbols that correspond to more than 1 gene_id. 
#-As the rest of the processing deals focuses on merging via gene_id, and not symbol
#-For example: Aldoa,
#http://uswest.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000114515;r=7:126399962-126408280
#http://uswest.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000030695;r=7:126394406-126399923
#-Both are the same gene symbol, but have different gene IDs. 
#-For simplicity sake, we will remove Gene_IDs until each Symbol has only 1 Gene_ID

# these 20 Symbols were left as is associated with their various unique gene ids.
pc_sub_duplicates <- pc %>% 
  select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  filter(Symbol !="" & duplicated(Symbol)) %>%
  distinct(Symbol)

#---------------------------------------------------------------------------
#Merge pc table onto the expression test frame
#---------------------------------------------------------------------------
# Protein coding (pc) table must be merged onto the test frame so that
# protein coding genes can be identified from our measured RNAseq counts.

test_merged <- test_table %>%
  left_join(y = pc_sub, by = "Gene_ID") %>% 
  relocate(Symbol, .after = Gene_ID)

#---------------------------------------------------------------------------
#Looking at breakdown of protein coding genes that did not join on to the count table
#---------------------------------------------------------------------------
# A quick glance at the protein coding genes that were not identified 
# in the count table shows that many were
# Olfactory proteins, and RIKEN-named genes
did_not_join <- filter(pc_sub, !pc_sub$Gene_ID %in% test_merged$Gene_ID)

#---------------------------------------------------------------------------
#Process the merged data frame. Into gene_ids with symbols(pc), and gene_ids without symbols
#---------------------------------------------------------------------------
#We want to identify the protein coding genes that are well documented;
#thus we are filtering for genes with symbols.
#We also will identify the Gene_IDs without symbols just to identify and explore

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

#---Save with_Symbols as a TSV

#test_table_ID is the test table's actual sample ID. Using the ID in the
#saved products maks it clearer what is being saved
test_table_ID <- str_replace(names(l_expr_tables[index]), pattern = ".tsv", replacement = "")
write_tsv(with_symbols, file = paste0("Data/test_table_analysis/",test_table_ID,"_pc_counts.tsv"))

#---------------------------------------------------------------------------
# Generating Alex-requested summary stats
#---------------------------------------------------------------------------

n_matched_to_symbol <- nrow(with_symbols)

frac_matched_to_symbol <- nrow(with_symbols)/nrow(test_merged) 

n_not_matched_to_symbol <- nrow(empty_symbols)

frac_not_matched_to_symbol <-  nrow(empty_symbols)/nrow(test_merged)

stopifnot(frac_matched_to_symbol+frac_not_matched_to_symbol == 1)

n_not_matched_pc_table <- nrow(pc_sub) - sum(pc_sub$Gene_ID %in% test_merged$Gene_ID)

#pc_sub has no empty rows. Each match has a Symbol
#rows removed in n_not_matched_to_symbol did not have a symbol. 
# Likely non-coding genes and such

#---------------------------------------------------------------------------
# This part of the script is for plotting
#---------------------------------------------------------------------------
names <- factor(x = c("Size of PC Table", 
                      "Num Gene Id's that didn't match", 
                      "Num Gene Id's that matched", 
                      "Size of Merged Count Table", 
                      "Num of Measured Proteins with Symbols", 
                      "Num of Measured Proteins without Symbols"),
                levels = c("Size of PC Table", 
                           "Num Gene Id's that didn't match", 
                           "Num Gene Id's that matched", 
                           "Size of Merged Count Table", 
                           "Num of Measured Proteins with Symbols", 
                           "Num of Measured Proteins without Symbols"))

lengths <- c(nrow(pc_sub), 
             n_not_matched_pc_table, 
             sum(pc_sub$Gene_ID %in% test_merged$Gene_ID),
             nrow(test_merged), 
             n_matched_to_symbol, 
             n_not_matched_to_symbol)

df_summary_stats <- data.frame(x = names, y = lengths)


pl_summary_stats <- ggplot(df_summary_stats) +
  geom_col(mapping = aes(x = names, y = lengths))+
  labs(title = test_table_ID)+
  xlab("Summary Stats")+
  ylab("Counts of Gene Id\'s")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
pl_summary_stats

write_delim(df_summary_stats,
            file = paste0("Data/test_table_analysis/", test_table_ID,"_df.tsv"),
            delim = "\t" )
ggsave(file = paste0("Data/test_table_analysis/"
                     ,test_table_ID,"_plot.png"),
       plot = pl_summary_stats)


#---------------------------------------------------------------------------
########################## 3) Generate  sample - expression data matrix
#--------------------------------------------------------------------------

# This script processes ALL count tables within l_expr_tables. 
# remove gene id versions from test data frames
# merge test data frame with pc table.
# Data frames are used to generate a  sample - expression data matrix
# sample - expression data matrix is saved


#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

# l_expr_tables is already opened: It is a list containing all of the count tables
# pc_sub is already generated: It is a processed pc table

#---Open meta_data
meta_data <- read.delim("Data/complete_meta_data.tsv",
                        stringsAsFactors = FALSE,
                        sep = "\t")

#---------------------------------------------------------------------------
# Prep PC and l_expr_tables for merging with each other, then merge
#--------------------------------------------------------------------------

#---Run remove_gene_id_vers on count tables to remove gene id versions
stopifnot("Gene_ID" %in% colnames(l_expr_tables[[1]]))
l_expr_tables_noid <- lapply(l_expr_tables, remove_gene_id_vers)

#---merge pc_sub table with all count tables
l_expr_merged <- lapply(l_expr_tables_noid, merge_count_pc, pc_sub)

#---------------------------------------------------------------------------
# Get only Protein Coding genes from merged expression - pc data frames. 
#---------------------------------------------------------------------------
# Now that our count tables have Protein Symbols, we can 
# Use the presence of a symbol to filter for protein coding genes 
l_expr_symbols <- lapply(l_expr_merged, get_genes_with_symbols)

#---------------------------------------------------------------------------
# Creating the sample - expression data matrix
#--------------------------------------------------------------------------
#---Test to see that each frame in l_expr_symbols has the same order of symbols
# If this is true, then we can just iteratively add the lists onto a matrix

stopifnot(test_orders(l_expr_symbols))


#---Create empty expression matrix and fill with expression data
count_matrix <- matrix(data = 0,
                            nrow = nrow(l_expr_symbols$ENCFF227HKF),
                            ncol = length(l_expr_symbols))

rownames(count_matrix) <- l_expr_symbols$ENCFF227HKF$Symbol
colnames(count_matrix) <- names(l_expr_symbols)

for (colname in colnames(count_matrix)) { 
  count_matrix[ , colname] <- l_expr_symbols[[colname]][["TPM"]]
}

#---Save expression matrix as rds

saveRDS(object = count_matrix, file = "Data/pc_count_matrix.rds")

#---------------------------------------------------------------------------
# Building an AVG Matrix
#---------------------------------------------------------------------------
# As seen in meta_data, every mouse developmental stage has 2 samples assosiated
# with it. We want to create a new matrix with those replicates averaged so
# we can easily perform analysis on how our counts change over the dev stages
# Build a new matrix containing the averages of counts across replicates
# Called avg_count_matrix

uniq_dev_stages <- unique(meta_data$dev_stage)

avg_count_matrix <- matrix(data = 0,
                           nrow = nrow(count_matrix),
                           ncol = ncol(count_matrix)/2)
rownames(avg_count_matrix) <- rownames(count_matrix)
colnames(avg_count_matrix) <- uniq_dev_stages 

for (mydev_stage in uniq_dev_stages) { 
  meta_subset <- filter(meta_data, meta_data$dev_stage == mydev_stage)
  avg <- rowMeans(count_matrix[,meta_subset$id])
  avg_count_matrix[,mydev_stage] <- avg
}

#---Check to confirm the averaging functioned properly

for (check_number in 1:5) {
  check_name <- sample(rownames(avg_count_matrix), size = 1)
  
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,1:2]),
                                  avg_count_matrix[check_name, 1]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,3:4]),
                                  avg_count_matrix[check_name, 2]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,9:10]),
                                  avg_count_matrix[check_name, 5]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,15:16]),
                                  avg_count_matrix[check_name, 8]))
}

#---Save avg expression matrix as rds

saveRDS(object = avg_count_matrix, file = "Data/avg_pc_count_matrix.rds")


#---------------------------------------------------------------------------
# Creating a difference count matrix
#---------------------------------------------------------------------------
# We want to identify which genes in a pair of replicates differ the most
# in count value.


uniq_dev_stages <- unique(meta_data$dev_stage)

diff_count_matrix <- matrix(0, 
                            nrow = nrow(avg_count_matrix),
                            ncol = ncol(avg_count_matrix))
rownames(diff_count_matrix) <- rownames(avg_count_matrix)
colnames(diff_count_matrix) <- colnames(avg_count_matrix)

for (mydev_stage in uniq_dev_stages) { 
  meta_subset <- filter(meta_data, meta_data$dev_stage == mydev_stage)
  stopifnot(length(meta_subset$id) == 2)
  diff_count_matrix[, mydev_stage] <-  abs(count_matrix[ , meta_subset$id[1]] - count_matrix[ , meta_subset$id[2]])
}


#---Check to confirm the averaging functioned properly

for (check_number in 1:5) {
  check_name <- sample(rownames(diff_count_matrix), size = 1)
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 1] - count_matrix[check_name, 2])
                                  , diff_count_matrix[check_name, 1]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 3] - count_matrix[check_name, 4])
                                  , diff_count_matrix[check_name, 2]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 9] - count_matrix[check_name, 10])
                                  , diff_count_matrix[check_name, 5]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 15] - count_matrix[check_name, 16])
                                  , diff_count_matrix[check_name, 8]))
  
  stopifnot(assertthat::are_equal(abs(count_matrix[check_name, 13] - count_matrix[check_name, 14])
                                  , diff_count_matrix[check_name, 7]))
  
}


#---Save avg expression matrix as rds

saveRDS(object = diff_count_matrix, file = "Data/diff_pc_count_matrix.rds")



