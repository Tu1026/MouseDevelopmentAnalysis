# This script processes 1 merged expression - pc data frame. It is named as test_merged as an RDS object.
#

source(file = "Main/functions.R")

library(skimr)
library(tidyverse)

#---Open up test_merged rds object

test_merged <- readRDS(file = "Data/test_merged")

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
tryCatch(
  {
    pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
  },
  error = function(cond) {
    message(cond)
    message("Couldn't open the pc table. It should be in Data/")
  }
)

pc_length <- nrow(pc)

tryCatch(
  {
    test_frame <- readRDS("Data/test_frame")
  },
  error = function(cond) {
    message(cond)
    message("Couldn't open test_frame. It should be in Data/")
  }
)

test_frame_length <- nrow(test_frame)

names <- c("Protein Coding Table", "Expression Table", "Merged Frame", "Transcripts with Symbol",
           "Transcripts with no Symbol", "Unique Expressed Genes")
names <- factor(names, levels = names)

lengths <- c(nrow(pc), nrow(test_frame), nrow(test_merged), n_pc_expr_1_Sym, n_pc_expr_1_noSym, nrow(test_merged_uniq))

summary_frame <- data_frame("names" = names, 
                            "lengths" = lengths)
summary_plot <- ggplot(summary_frame, mapping = aes(x = names, y = lengths)) +
  geom_col()

#---Save plot as rds object

if (!dir.exists("Plots")) {
  dir.create("Plots")
}

saveRDS(summary_plot, file = "Plots/summary_plot")




#---Archived plotting
# #---Create Summary Data Frame for plotting
# df_summary <- data.frame(
#   "Merged Transcripts" = nrow(pc_expr_1),
#   "Transcripts with Symbol(PC Genes)" =n_pc_expr_1_Sym,
#   "Transcripts with no Symbol" = n_pc_expr_1_noSym,
#   "Expressed PC Transcripts" = n_count_1_expr_1_Sym,
#   "Non Expressed PC Transcripts" = n_count_0_expr_1_Sym,
#   "Expressed Genes" = count_1_pc_genes,
#   "Non Expressed Genes" = count_0_pc_genes
# )
# 
# row1 <- c("Merged Transcripts","Transcripts with Symbol(PC Genes)","Transcripts with no Symbol" , "Expressed PC Transcripts", "Non Expressed PC Transcripts",
#           "Expressed Genes", "Non Expressed Genes")
# row2 <- c(nrow(pc_expr_1), n_pc_expr_1_Sym, n_pc_expr_1_noSym, n_count_1_expr_1_Sym, n_count_0_expr_1_Sym, count_1_pc_genes, count_0_pc_genes)
# df_melted <- data.frame(names = factor(row1, levels = row1), 
#                         values = row2)
# #TODO melt this data frame properly
# 
# summary_plot <- ggplot(df_melted) + 
#   geom_col(mapping = aes(x= names, y = values))
# 
# summary_plot



#---Archive

# #---Filter For the protein coding transcripts with expected count 0
# count_0_expr_1_Sym <- pc_expr_1_Sym %>%
#   filter(pc_expr_1_Sym$expected_count == 0)
# 
# #---Filter For the number of protein coding transcripts with expected count 0
# n_count_0_expr_1_Sym <- count_0_expr_1_Sym %>%
#   nrow()
# 
# #---Filter for Protein coding transcripts that have expected count >0
# count_1_expr_1_Sym <- pc_expr_1_Sym %>%
#   filter(pc_expr_1_Sym$expected_count>0)
# 
# #---Filter for the number of Protein coding transcripts that have expected count >0
# n_count_1_expr_1_Sym <- count_1_expr_1_Sym %>%
#   nrow()
# 
# #---Check to make sure the number of transcripts with expression 0 + number of transcripts with expression >0 is equal to the total number of PC transcripts
# stopifnot(n_count_1_expr_1_Sym+n_count_0_expr_1_Sym == n_pc_expr_1_Sym)
# 
# #---Filter for the number of unique pc transcripts in total (expressed + non expressed)
# pc_expr_1_Sym_uniq <- pc_expr_1_Sym%>%
#   group_by(Symbol) %>%
#   tally()%>%
#   nrow()
# 
# #---For Expressed PC Transcripts. Get the number of Symbols (Unique Genes)
# count_1_pc_genes <- count_1_expr_1_Sym %>%
#   group_by(Symbol) %>%
#   tally()%>%
#   nrow()
# 
# #---Filter for the number of unique genes NOT expressed
# count_0_pc_genes <- count_0_expr_1_Sym %>%
#   group_by(Symbol) %>%
#   tally()%>%
#   nrow()
# 
# #---Check to make sure the number of unique genes expressed + number of unique genes not expressed == number of unique transcripts with Symbols
# stopifnot(count_0_pc_genes+count_1_pc_genes==pc_expr_1_Sym_uniq)
# 
# #---Rename count_1_expr_1_Sym to "Expressed Genes" Because thats what it is
# expressed_genes <- count_1_expr_1_Sym %>%
#   group_by(Symbol)

