source(file = "Main/functions.R")
source("Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")
source ("Main/05_merge.R")
source("Main/06_process_merged_tables.R")
source("Main/07_selecting_longest_transcript.R")

glimpse(expressed_genes)
glimpse(expressed_genes_pro_transcript)

#Lets inspect this data frame. Lets get some summary statistics from our df with uniquely expressed genes (no transcripts)

#---basic summary stats for expression_count
expression_summary_stats <- expressed_genes_pro_transcript %>%
  ungroup()%>%
  select(Symbol, expected_count)%>%
  skim()

#---Highest Expressed Genes using Expression_count
highest_expressions <- expressed_genes_pro_transcript %>%
  ungroup()%>%
  arrange(desc(expressed_genes_pro_transcript$expected_count))%>%
  select(Symbol, expected_count)%>%
  head(n = 100)

#---Lowest Expressed Genes using Expression count
lowest_expressions <- expressed_genes_pro_transcript %>%
  ungroup()%>%
  arrange(expressed_genes_pro_transcript$expected_count)%>%
  select(Symbol, expected_count)%>%
  head(n = 100)

#---Visualizing expression statistics using expression_count
expression_hist <- ggplot(data = expressed_genes_pro_transcript)+
  geom_histogram(mapping = aes(expressed_genes_pro_transcript$expected_count), bins = 100)+
  scale_x_log10() +
  labs(title = "Gene Expression") +
  xlab("Gene Expression Count")+
  ylab("n")
expression_hist


