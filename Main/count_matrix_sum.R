library(tidyverse)
library(ggplot2)

source('Main/functions.R')

#---------------------------------------------------------------------------
# Open Correlation Matrix and meta_data
#---------------------------------------------------------------------------

count_matrix <- readRDS(file = "Data/pc_count_matrix.rds")
avg_count_matrix <- readRDS(file = "Data/avg_pc_count_matrix.rds")
meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")

#---------------------------------------------------------------------------
# Removing 0s
#---------------------------------------------------------------------------
count_matrix <- count_matrix[ count_matrix[ , 1] > 0 |
                                count_matrix[ , 2] > 0 | 
                                count_matrix[ , 3] > 0 |
                                count_matrix[ , 4] > 0 |
                                count_matrix[ , 5] > 0 |
                                count_matrix[ , 6] > 0 |
                                count_matrix[ , 7] > 0 |
                                count_matrix[ , 8] > 0 |
                                count_matrix[ , 9] > 0 | 
                                count_matrix[ , 10] > 0 |
                                count_matrix[ , 11] > 0 |
                                count_matrix[ , 12] > 0 |
                                count_matrix[ , 13] > 0 |
                                count_matrix[ , 14] > 0 |
                                count_matrix[ , 15] > 0 |
                                count_matrix[ , 16] > 0 , ]

avg_count_matrix <- avg_count_matrix[ avg_count_matrix[ , 1] > 0 |
                                        avg_count_matrix[ , 2] > 0 | 
                                        avg_count_matrix[ , 3] > 0 |
                                        avg_count_matrix[ , 4] > 0 |
                                        avg_count_matrix[ , 5] > 0 |
                                        avg_count_matrix[ , 6] > 0 |
                                        avg_count_matrix[ , 7] > 0 |
                                        avg_count_matrix[ , 8] > 0 , ]

#---------------------------------------------------------------------------
# Log2 normalizing 
#---------------------------------------------------------------------------
count_matrix <- log2(count_matrix+1)
avg_count_matrix <- log2(avg_count_matrix+1)

#---------------------------------------------------------------------------
# Opening Expr Tables
#---------------------------------------------------------------------------

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
#--------------------------------------------------------------------------
target_directory <- "Data/Count_tables"
l_expr_tables <- lapply(ordered_meta_data_id_names,
                        open_expr_table,
                        target_directory)

names(l_expr_tables) <- str_replace(ordered_meta_data_id_names, ".tsv","")



#---------------------------------------------------------------------------
# Graphing
#---------------------------------------------------------------------------


df_summary_stats <- data.frame( names = c("Annotated Genes", "Measured PC Genes"),
                                lengths = c(nrow(l_expr_tables[[1]]), nrow(count_matrix)))

pl_summary_stats <- ggplot(df_summary_stats) +
  geom_col(mapping = aes(x = names, y = lengths), fill = "steelblue")+
  labs(title = "Measured Genes")+
  xlab("Summary Stats")+
  ylab("Counts")+
  theme(axis.text.x = element_text(angle = 0))+
  theme_classic()
pl_summary_stats




#---------------------------------------------------------------------------
# save count_matrix with 0s
#---------------------------------------------------------------------------

saveRDS(count_matrix, file = "Data/pc_count_matrix0slog2.rds")
saveRDS(avg_count_matrix, file = "Data/avg_pc_count_matrix0slog2.rds")


#---------------------------------------------------------------------------
# count_matrix histogram that doenst have 0s removed
#---------------------------------------------------------------------------

count_matrix_all <- readRDS(file = "Data/pc_count_matrix.rds")
count_matrix_all <- log2(count_matrix_all+1)

df_for_plot <- as.data.frame(rowMeans(count_matrix_all))
colnames(df_for_plot) <- "Mean Expression"
view(df_for_plot)

ggplot(data = df_for_plot) + 
  geom_histogram( mapping = aes(df_for_plot[, 1]), fill = "steelblue") +
  labs(x = "Expression", y = "Count", title = "Mean Expression: Non-Expressed Genes Retained")

#---------------------------------------------------------------------------
# count_matrix histogram that has 0s removed
#---------------------------------------------------------------------------

df_for_plot <- as.data.frame(rowMeans(count_matrix))

ggplot(data = df_for_plot) + 
  geom_histogram( mapping = aes(df_for_plot[,1]), fill = "steelblue") +
  labs(x = "Expression", y = "Count", title = "Mean Expression: Non-Expressed Genes Removed")


