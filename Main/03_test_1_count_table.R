

#---------------------------------------------------------------------------
########################## Processing 1 count table
#--------------------------------------------------------------------------

# Process a single data frame within l_expr_tables. 
# Remove gene id versions from test data frame
# Merge test data frame with pc table. Name it test_merged
# It then processes the test data frame. This processing will create a 
# df containing Measured genes
# Get some visualization about the filtering is generated.
# Get some summary statistics about our data frame


source('Main/functions.R')

library(tidyverse)
library(assertthat)
library(stringr)
library(skimr)

#---------------------------------------------------------------------------
#Setup and opening files
#---------------------------------------------------------------------------


#---Choose which count table to process
# index is the frame that you are selecting from l_expr_tables
count_tables_names <- list.files("Data/Count_tables")
index <- 1
test_table <- open_expr_table(target_directory = "Data/Count_tables",
                              expression_table.tsv = count_tables_names[index])

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
test_table_ID <- str_replace(count_tables_names[index], pattern = ".tsv", replacement = "")
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
names <- factor(x = c(
                      "Annotated Genes", 
                      "Measured PC Genes"),

                levels = c(
                           "Annotated Genes", 
                           "Measured PC Genes")
                )
lengths <- c( 
             nrow(test_table),
             nrow(count_matrix))

df_summary_stats <- data.frame(x = names, y = lengths)

### Open meta_data for the naminf of this plot
fb_meta <- read.delim(file = "Data/complete_meta_data.tsv",
           sep = "\t",
           stringsAsFactors = FALSE)
test_table_stage <- fb_meta %>%
  filter(id == test_table_ID) %>%
  select(dev_stage)

pl_summary_stats <- ggplot(df_summary_stats) +
  geom_col(mapping = aes(x = names, y = lengths), fill = "steelblue")+
  labs(title = paste0(test_table_ID, " (",test_table_stage$dev_stage, ") "))+
  xlab("Summary Stats")+
  ylab("Counts of Gene Id\'s")+
  theme(axis.text.x = element_text(angle = 0))+
  theme_classic()
pl_summary_stats

write_delim(df_summary_stats,
            file = paste0("Data/test_table_analysis/", test_table_ID,"_df.tsv"),
            delim = "\t" )
ggsave(file = paste0("Data/test_table_analysis/"
                     ,test_table_ID,"_plot.png"),
       plot = pl_summary_stats,
       dpi = 1000,
       width = 5.5, 
       height = 4.5)

