library(tidyverse)
library(assertthat)
library(stringr)
library(Hmisc)
library(corrplot)
# library(preprocessCore)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

source('Main/functions.R')


#---------------------------------------------------------------------------
# Open Correlation Matrix and meta_data
#---------------------------------------------------------------------------

count_matrix <- readRDS(file = "Data/pc_count_matrix.rds")
avg_count_matrix <- readRDS(file = "Data/avg_pc_count_matrix.rds")
diff_count_matrix <- readRDS(file = "Data/diff_pc_count_matrix.rds")
meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")


#---------------------------------------------------------------------------
# createdir
#---------------------------------------------------------------------------
dir.create(path = "Data/PCA", showWarnings = FALSE)

#---------------------------------------------------------------------------
# PCA
#---------------------------------------------------------------------------

#prcomp expects genes to be columns and samples to be rwos so we have to t()
pca <- prcomp(t(count_matrix))

# x contains the principle components ( the vecrtors of length1)
plot(pca$x[,1], pca$x[,2])

#there should be 16 principle components, 1 for each sample
plot(pca$x[,1], pca$x[,3])

plot(pca$x[,2], pca$x[,3])

# To identify how much of variation in the data a PC can account for
# we square sdev

pca_var <- pca$sdev^2

# To get the percentage value we have to devide by total variation
pca_var_per <- round(pca_var/sum(pca_var) *100, 1)
pca_var_per[1]

#to get a scree plot we just have to plot this

barplot(pca_var_per)

# GGplot plotting
dev_stages <- c(rep("E10.5",2) ,
               rep("E11.5",2) ,
               rep("E12.5",2) ,
               rep("E13.5",2) ,
               rep("E14.5",2) ,
               rep("E15.5",2) ,
               rep("E16.5",2) ,
               rep("P0",2)
               )

df_pca_val <- data.frame(Dev_Stages = dev_stages, 
                         pca1 = pca$x[,1],
                         pca2 = pca$x[,2]
                         )
pca_plot <- ggplot(df_pca_val) +
  geom_point(mapping = aes(x=pca1,
                           y=pca2,
                           color = Dev_Stages)) + 
  ggtitle("PCA plot") + 
  xlab("PCA1")+
  ylab("PCA2")

pca_plot

ggsave(filename = "Data/PCA/pca_plot.png")

df_pca_var <- data.frame(names = factor(colnames(pca$x),
                                        levels = colnames(pca$x)),
                     per_var <- pca_var_per)
                     
                     
scree <- ggplot(df_pca_var) + 
  geom_col(mapping = aes(x = names,
                         y = per_var), fill = "steelblue")+
  ggtitle("Scree plot")+
  xlab("PCA component")+
  ylab("% Variation accounted for")
  

scree

ggsave(filename = "Data/PCA/scree.png")


