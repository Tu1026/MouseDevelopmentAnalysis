source(file = "Main/functions.R")
source("Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")
source ("Main/05_merge.R")
source("Main/06_process_merged_tables.R")

#TODO. I'm confused here. I thought we were taking the transcript that is the longest for a given gene

#This script is for determining what transcript ID we want to use as our gene

#---Our Genes are the ones with the highest expression #---Note: They all have the same expressiion because of the merge
expressed_genes_pro_max <- expressed_genes %>%
  group_by(Symbol) %>%
  top_n(n=1)


#---Selecting the genes with the longest Transcript
expressed_genes_pro_transcript <- expressed_genes %>%
  ungroup() %>%
  mutate("len_Transcript_ID" = expressed_genes$End - expressed_genes$Start)%>%
  group_by(Symbol)%>%
  top_n(n=1)
