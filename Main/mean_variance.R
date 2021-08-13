library(tidyverse)
library(assertthat)
library(stringr)
library(Hmisc)
library(corrplot)
# library(preprocessCore)
library(cowplot)
library(pheatmap)
library(RColorBrewer)

source('Main/functions.R')


#---------------------------------------------------------------------------
# Open Correlation Matrix and meta_data
#---------------------------------------------------------------------------

count_matrix <- readRDS(file = "Data/pc_count_matrix.rds")
avg_count_matrix <- readRDS(file = "Data/avg_pc_count_matrix.rds")
diff_count_matrix <- readRDS(file = "Data/diff_pc_count_matrix.rds")
meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")

#---------------------------------------------------------------------------
# Log 10+1 normalizing 
#---------------------------------------------------------------------------
#no normalization here

#---------------------------------------------------------------------------
# Create place for new matrix exploration graphs 
#---------------------------------------------------------------------------
dir.create(path = "Data/Mean_Variance", showWarnings = FALSE)

#---------------------------------------------------------------------------
# Create place for new matrix exploration graphs 
#---------------------------------------------------------------------------
row_var <- rowVar(count_matrix)

row_sum <- rowSums(count_matrix)
row_mean <- rowMeans(count_matrix)
row_stdev <- rowSD(count_matrix)
row_rel_sd <- row_stdev/row_mean


summary_matrix <- matrix(c(row_sum,
                           row_mean,
                           row_stdev,
                           row_rel_sd),
                         nrow = length(row_sum),
                         ncol = 4)
colnames(summary_matrix) <-  c("row_sum",
                               "row_mean",
                               "row_stdev",
                               "row_rel_sd")
rownames(summary_matrix) <-  names(row_sum)

summary_matrix

#---------------------------------------------------------------------------
# mean variance plot for genes
#---------------------------------------------------------------------------

plot(row_mean, sqrt(row_var))

stopifnot(FALSE)

row1 <- count_matrix[1,]
var(x = row1)
?var
