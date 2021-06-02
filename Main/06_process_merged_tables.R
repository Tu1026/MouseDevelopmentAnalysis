source(file = "Main/functions.R")
source("Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")
source ("Main/05_merge.R")

library(skimr)

#In this script I'm going to try and process 1 single l_pc_expr data frame
#The end result is a data frame with summary information as well as some basic visualization 

pc_expr_1 <- l_pc_expr[[1]]

#---Filter out transcripts that do have a protein symbol

pc_expr_1_Sym <- pc_expr_1 %>%
  filter(!(is.na(pc_expr_1$Symbol)))

#---Get the number of transcripts with symbols
n_pc_expr_1_Sym <- pc_expr_1_Sym %>%
  nrow()

#---Check that all of your transcriopts do indeed have symbols
stopifnot(!(all(is.na(pc_expr_1_Sym$Symbol))))


#---Filter out transcripts that NOT have a protein symbol
pc_expr_1_noSym <- pc_expr_1 %>%
  filter((is.na(pc_expr_1$Symbol)))

#---Get the number of transcripts with symbols
n_pc_expr_1_noSym <- pc_expr_1_noSym %>%
  nrow()

#---Check that all of your transcriopts do  NOT have symbols
stopifnot(!(all(is.na(pc_expr_1_noSym$Symbol))))

#Check to make sure your number of symbol transcripts + non-symbol transcripts adds up to total
stopifnot(n_pc_expr_1_noSym+n_pc_expr_1_Sym == nrow(pc_expr_1))

#---Filter For the protein coding transcripts with expected count 0
count_0_expr_1_Sym <- pc_expr_1_Sym %>%
  filter(pc_expr_1_Sym$expected_count == 0)

#---Filter For the number of protein coding transcripts with expected count 0
n_count_0_expr_1_Sym <- count_0_expr_1_Sym %>%
  nrow()

#---Filter for Protein coding transcripts that have expected count >0
count_1_expr_1_Sym <- pc_expr_1_Sym %>%
  filter(pc_expr_1_Sym$expected_count>0)

#---Filter for the number of Protein coding transcripts that have expected count >0
n_count_1_expr_1_Sym <- count_1_expr_1_Sym %>%
  nrow()

#---Check to make sure the number of transcripts with expression 0 + number of transcripts with expression >0 is equal to the total number of PC transcripts
stopifnot(n_count_1_expr_1_Sym+n_count_0_expr_1_Sym == n_pc_expr_1_Sym)

#---Filter for the number of unique pc transcripts in total (expressed + non expressed)
pc_expr_1_Sym_uniq <- pc_expr_1_Sym%>%
  group_by(Symbol) %>%
  tally()%>%
  nrow()

#---For Expressed PC Transcripts. Get the number of Symbols (Unique Genes)
count_1_pc_genes <- count_1_expr_1_Sym %>%
  group_by(Symbol) %>%
  tally()%>%
  nrow()

#---Filter for the number of unique genes NOT expressed
count_0_pc_genes <- count_0_expr_1_Sym %>%
  group_by(Symbol) %>%
  tally()%>%
  nrow()

#---Check to make sure the number of unique genes expressed + number of unique genes not expressed == number of unique transcripts with Symbols
stopifnot(count_0_pc_genes+count_1_pc_genes==pc_expr_1_Sym_uniq)

#---Rename count_1_expr_1_Sym to "Expressed Genes" Because thats what it is
expressed_genes <- count_1_expr_1_Sym %>%
  group_by(Symbol)



#---Create Summary Data Frame for plotting
df_summary <- data.frame(
  "Merged Transcripts" = nrow(pc_expr_1),
  "Transcripts with Symbol(PC Genes)" =n_pc_expr_1_Sym,
  "Transcripts with no Symbol" = n_pc_expr_1_noSym,
  "Expressed PC Transcripts" = n_count_1_expr_1_Sym,
  "Non Expressed PC Transcripts" = n_count_0_expr_1_Sym,
  "Expressed Genes" = count_1_pc_genes,
  "Non Expressed Genes" = count_0_pc_genes
)

row1 <- c("Merged Transcripts","Transcripts with Symbol(PC Genes)","Transcripts with no Symbol" , "Expressed PC Transcripts", "Non Expressed PC Transcripts",
          "Expressed Genes", "Non Expressed Genes")
row2 <- c(nrow(pc_expr_1), n_pc_expr_1_Sym, n_pc_expr_1_noSym, n_count_1_expr_1_Sym, n_count_0_expr_1_Sym, count_1_pc_genes, count_0_pc_genes)
df_melted <- data.frame(names = factor(row1, levels = row1), 
                        values = row2)
#TODO melt this data frame properly

summary_plot <- ggplot(df_melted) + 
  geom_col(mapping = aes(x= names, y = values))

summary_plot

