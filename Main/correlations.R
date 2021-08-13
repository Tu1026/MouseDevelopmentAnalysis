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
# Generating Pearson Correlation Matrix
#---------------------------------------------------------------------------
# Generating pearson correlation matrix to see correlation between sample's gene
# counts

# Also generate correlation matrix's p-values

pear_rcorr_object <- rcorr(count_matrix, type = "pearson")
pear_matrix <- pear_rcorr_object$r
pear_matrix_pvalues <- pear_rcorr_object$P

#---Generating heatmap for sample correlation matrix

#since the heatmap is mirrored, we will make it a triangle heatmap

# diag(pear_matrix) <- NA # remove the diagnomal row of the matrix
pear_matrix[(upper.tri(pear_matrix))] <- NA # Remove the entire upper triangle

pear_heatmap <- pheatmap(pear_matrix, 
                         main = "Pearson Sample Correlation Matrix",
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         breaks = rev(c(1.0,0.99999999, 0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6)),
                         color = rev(brewer.pal(n = 10, name ="RdYlBu")),
                         legend_breaks = c(1.0,0.9,0.8,0.7,0.6),
                         na_col = "white",
                         border_color = "white",
                         labels_row = c("E10.5-1",
                                        "E10.5-2",
                                        "E11.5-1",
                                        "E11.5-2",
                                        "E12.5-1",
                                        "E12.5-2",
                                        "E13.5-1",
                                        "E13.5-2",
                                        "E14.5-1",
                                        "E14.5-2",
                                        "E15.5-1",
                                        "E15.5-2",
                                        "E16.5-1",
                                        "E16.5-2",
                                        "P0-1",
                                        "P0-2"),
                         labels_col = c("E10.5-1",
                                        "E10.5-2",
                                        "E11.5-1",
                                        "E11.5-2",
                                        "E12.5-1",
                                        "E12.5-2",
                                        "E13.5-1",
                                        "E13.5-2",
                                        "E14.5-1",
                                        "E14.5-2",
                                        "E15.5-1",
                                        "E15.5-2",
                                        "E16.5-1",
                                        "E16.5-2",
                                        "P0-1",
                                        "P0-2")
)

# --- Saving 

dir.create(path = "Data/correlations")

# 
# ggsave(file = paste0("Data/correlations/sample_pearson"),
#        plot = pear_heatmap)

pear_heatmap


#---------------------------------------------------------------------------
# Generating spearman Correlation Matrix
#---------------------------------------------------------------------------
# Generating spearman correlation matrix to see correlation between sample's gene
# counts if they were ordered

# Also generate correlation matrix's p-values

spear_rcorr_object <- rcorr(count_matrix, type = "spearman")
spear_matrix <- spear_rcorr_object$r
spear_matrix_pvalues <- spear_rcorr_object$P

#---Generating heatmap for sample correlation matrix

#since the heatmap is mirrored, we will make it a triangle heatmap

spear_matrix[(upper.tri(spear_matrix))] <- NA # Remove the entire upper triangle

spear_heatmap <- pheatmap(spear_matrix, 
                         main = "Spearman Sample Correlation Matrix",
                         cluster_cols = FALSE,
                         cluster_rows = FALSE,
                         breaks = rev(c(1.0,0.99999999, 0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6)),
                         color = rev(brewer.pal(n = 10, name ="RdYlBu")),
                         legend_breaks = c(1.0,0.9,0.8,0.7,0.6),
                         na_col = "white",
                         border_color = "white"
                         
)


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

#---------------------------------------------------------------------------
# Comments on Spearman and Pearson Correlation Matrices
#---------------------------------------------------------------------------

# The pearson correlation decreases as the mouse samples progress through
# the dev stages

# On the other hand, the spearman correlation stays relatively constant
# as the mouse samples progress through the dev stages

# Whereas the pearson is correlating pair-wise values which are paired with 
# respect to the Gene Symbols. It also assumes that both samples are normally
# distributed. The decrease of correlation coefficients across dev stages
# suggests that expression is changing. However the high correlation between
# replicates and close dev stages suggests the quality of the replicates
# and data is strong

# The high spearmen correlation coefficients suggest that all samples
# have a storng monotonic relationship. 


#---------------------------------------------------------------------------
# Comparing the count values of two replicate samples
#---------------------------------------------------------------------------

# Generate a scatterplot that compares the counts of two replicate samples

Replicate_1 <- count_matrix[,1]
Replicate_2 <- count_matrix[,2]

replicate_df <- data.frame("ENCFF302TQO (E10.5-1)" = Replicate_1,
                           "ENCFF416BXV (E10.5-2)" = Replicate_2, 
                           row.names = rownames(Replicate_1))

replicate_scatterplot <- ggplot(replicate_df)+ 
  geom_point(mapping = aes(x = Replicate_1, y = Replicate_2))+
  labs( title = "Pair-wise counts of two replicates")+ 
  xlab("ENCFF302TQO (E10.5-1)")+
  ylab("ENCFF416BXV (E10.5-2)")

replicate_scatterplot

#---------------------------------------------------------------------------
# Comparing the count values of two non-replicate samples
#---------------------------------------------------------------------------

# Generate a scatterplot that compares the counts of two replicate samples

non_Replicate_1 <- count_matrix[,1]
non_Replicate_2 <- count_matrix[,14]

non_replicate_df <- data.frame(non_Replicate_1 = non_Replicate_1,
                               non_Replicate_2 = non_Replicate_2, 
                               row.names = rownames(Replicate_1))

non_replicate_scatterplot <- ggplot(non_replicate_df)+ 
  geom_point(mapping = aes(x = non_Replicate_1, y = non_Replicate_2))+
  labs( title = "Pair-wise counts of two non replicate samples")+
  xlab("ENCFF302TQO (E10.5-1)")+
  ylab("ENCFF895JXR (E14.5-2)")
non_replicate_scatterplot
ggsave(filename = "Data/correlations/non_replicate_scatterplot.png")


#---------------------------------------------------------------------------
# Comments on replicate scatterplot and non-replicate scatterplot
#---------------------------------------------------------------------------

# As expected, the scatterplot of the replicates showed a strong correlation.
# Whereas the non-replicate scatterplot was much more scattered,
# indicating a weaker correlation
# This is indicative of the calculated pearson correlation coefficients,
# which were 0.99 and 0.60 for the two graphs respectively
