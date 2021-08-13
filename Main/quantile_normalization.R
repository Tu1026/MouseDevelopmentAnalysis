# install.packages("preprocessCore")

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
# Do the samples need to be quantile normalized
#---------------------------------------------------------------------------


plot(density(count_matrix[,1]),
     main = "Density of E10.5-1",
     xlab = "Density",
     ylab = "n")

plot(density(count_matrix[,1]),
     main = "Density of all Samples",
     xlab = "Density",
     ylab = "n")

for (i in 1:ncol(count_matrix)) {  # rest of the samples
  lines(density(count_matrix[, i]), col = "black")
}

# Density is very similar for all samples. Quantile normalization is not really 
# required

# For the rest of the exploration we will be working with the log values 
# of the count matrix