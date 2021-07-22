# install.packages("preprocessCore")

library(tidyverse)
library(assertthat)
library(stringr)
library(Hmisc)
library(corrplot)
# library(preprocessCore)
library(cowplot)

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

pear_heatmap <- heatmap(x = pear_matrix,
                   sym = TRUE,
                   Rowv = NA,
                   Colv = NA, 
                   revC= TRUE)

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

spear_heatmap <- heatmap(x = spear_matrix,
                        sym = TRUE,
                        Rowv = NA,
                        Colv = NA, 
                        revC= TRUE)


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

replicate_df <- data.frame(Replicate_1 = Replicate_1,
                           Replicate_2 = Replicate_2, 
                           row.names = rownames(Replicate_1))

replicate_scatterplot <- ggplot(replicate_df)+ 
  geom_point(mapping = aes(x = Replicate_1, y = Replicate_2))+
  labs( title = "Pair-wise counts of two replicate samples")

replicate_scatterplot


#---------------------------------------------------------------------------
# Comparing the count values of two non-replicate samples
#---------------------------------------------------------------------------

# Generate a scatterplot that compares the counts of two replicate samples

non_Replicate_1 <- count_matrix[,1]
non_Replicate_2 <- count_matrix[,16]

non_replicate_df <- data.frame(non_Replicate_1 = non_Replicate_1,
                               non_Replicate_2 = non_Replicate_2, 
                                row.names = rownames(Replicate_1))

non_replicate_scatterplot <- ggplot(non_replicate_df)+ 
  geom_point(mapping = aes(x = non_Replicate_1, y = non_Replicate_2))+
  labs( title = "Pair-wise counts of two non replicate samples")

non_replicate_scatterplot


#---------------------------------------------------------------------------
# Comments on replicate scatterplot and non-replicate scatterplot
#---------------------------------------------------------------------------

# As expected, the scatterplot of the replicates showed a strong correlation.
# Whereas the non-replicate scatterplot was much more scattered,
# indicating a weaker correlation
# This is indicative of the calculated pearson correlation coefficients,
# which were 0.99 and 0.60 for the two graphs respectively

#---------------------------------------------------------------------------
# Identifying differentially expressed genes between replicates
#---------------------------------------------------------------------------

# Ideally, replicates should have similar counts for all genes. 
# Here we will explore the genes that dont and add them to alist

# --- Looking at E105 first
E105 <- diff_count_matrix[, colnames(diff_count_matrix)[1]]

E105_df <- data.frame("E105" = E105)
E105_df <- arrange(E105_df, desc(E105_df$E105))

summary(E105_df)
summary(log10(E105_df+1))
E105_hist <- ggplot(log10(E105_df+1), aes(x = E105)) + 
  geom_histogram(bins=100) +
  labs( title = "E105 most differentially expresed genes")
E105_hist

l_most_diff_genes <- list("E105" = head(E105_df, n = 10))

# --- Doing the same for all dev stages
l_most_diff_genes <- list()


get_most_diff <- function(column) {
  my_df <- data.frame("Diff_count" = column)
  arranged_df <- arrange(my_df, desc(my_df))
  return(top_n(arranged_df, n = 10))
}

most_diff_genes <- apply(count_matrix, MARGIN = 2, FUN = get_most_diff)
most_diff_genes

#---------------------------------------------------------------------------
# Identifying most expressed genes across samples
#---------------------------------------------------------------------------

# Here we will identify the genes that are most expressed in each sample

# --- Looking at E105 first
E105 <- avg_count_matrix[, colnames(avg_count_matrix)[1]]

E105_df <- data.frame("E105" = E105)
E105_df <- arrange(E105_df, desc(E105_df$E105))

l_most_expr_genes <- list("E105" = head(E105_df, n = 10))

# --- Doing the same for all dev stages
l_most_expr_genes <- list()


get_most_expr<- function(column) {
  my_df <- data.frame("count" = column)
  arranged_df <- arrange(my_df, desc(my_df))
  return(top_n(arranged_df, n = 10))
}

most_expr_genes <- apply(count_matrix, MARGIN = 2, FUN = get_most_expr)


#---------------------------------------------------------------------------
# Comments on most expr and most diff genes
#---------------------------------------------------------------------------

# It looks like Tuba1a is the most expressed gene
# The problem is that many of the genes that are highly expressed
# are also seen to have the most deviation between replicates
# so we have to ask ourselves how reliable these counts are


#---------------------------------------------------------------------------
# Exploring genes that are never expressed
#---------------------------------------------------------------------------

# THere are some genes that never got counted across all dev stages
# Here I've arbitrarily chosen <=1 as never expressed. May chance later


never_expressed_genes <- avg_count_matrix[  avg_count_matrix[ , 1] <= 1 &
                        avg_count_matrix[ , 2] <= 1 & 
                        avg_count_matrix[ , 3] <= 1 &
                        avg_count_matrix[ , 4] <= 1 &
                        avg_count_matrix[ , 5] <= 1 &
                        avg_count_matrix[ , 6] <= 1 &
                        avg_count_matrix[ , 7] <= 1 &
                        avg_count_matrix[ , 8] <= 1 , ]



#---------------------------------------------------------------------------
# Do the samples need to be quantile normalized
#---------------------------------------------------------------------------

count_matrix <- log10(count_matrix+1)

for (i in 1:ncol(count_matrix)) {  # rest of the samples
  lines(density(count_matrix[, i]), col = "black")
}

# Density is very similar for all samples. Quantile normalization is not really 
# required

# For the rest of the exploration we will be working with the log values 
# of the count matrix

#---------------------------------------------------------------------------
# Summarizing Gene Expression across rows and constructing matrix
#---------------------------------------------------------------------------

summary <- summary(count_matrix[,])

row_sum <- rowSums(count_matrix)
summary(row_sum)

row_mean <- rowMeans(count_matrix)
summary(row_mean)

row_stdev <- rowSD(count_matrix)

row_rel_sd <- row_stdev/row_mean

summary_matrix <- matrix(c(row_sum, row_mean, row_stdev, row_rel_sd),
                         nrow = length(row_sum),
                         ncol = 4)
colnames(summary_matrix) <-  c("row_sum",
                               "row_mean",
                               "row_stdev",
                               "row_rel_sd")
rownames(summary_matrix) <-  names(row_sum)

#---------------------------------------------------------------------------
# mean variance plot for al lgenes
#---------------------------------------------------------------------------

plot(row_mean, row_stdev)

head(sort(row_rel_sd, decreasing = TRUE), n = 1000)

row_rel_sd
head(row_rel_sd)

row_stdev["Clec2g"]
row_mean["Clec2g"]
row_rel_sd["Clec2g"]

row_stdev["Clec2g"]/row_mean["Clec2g"]


# The genes with the highest relative standard deviation might be of the greatest
# interest as they might vary significantly. Indicating that they 
# change throughout development, which suggests that they are tightly regulated
# and important to development.

# Question: Why are they all converging at 4? ?????


#---------------------------------------------------------------------------
# Genes that are never expressed
#---------------------------------------------------------------------------

# If a gene's row summary is equal to 0, then we are saying 
# that it was never expressed


never_expressed_genes <- row_sum [row_sum == 0]

expressed_genes <- row_sum[row_sum != 0]
expressed_genes

expressed_row_mean <- row_mean [ row_sum !=0 ]
expressed_row_sd <- row_stdev [ row_sum !=0 ]
expressed_row_rel_sd <- row_rel_sd [ row_sum !=0]

#---------------------------------------------------------------------------
# mean variance plot for expressed genes
#---------------------------------------------------------------------------

plot(expressed_row_mean, expressed_row_sd)
#It actually doesn't change because I imagine the 0s are already taken out for
# previous mean variance plot anyways

df_rel_sd <- data.frame(expressed_row_rel_sd)
ggplot(df_rel_sd) + 
  geom_histogram(mapping = aes(expressed_row_rel_sd))


head(sort(expressed_row_rel_sd, decreasing = TRUE), n = 1000)

summary_matrix [names(sort(summary_matrix[,4], decreasing = TRUE)), ]

# Notably, all of the higehst mean variance genes have really low expression


#---------------------------------------------------------------------------
# Investigatging Clec2g
#---------------------------------------------------------------------------
# I've chosen this gene artibrarily. Just to see what it looks like over time

count_matrix["Clec2g",]

# as expected, there is very little expression over time. It does change
# but its hard to make any actual assumptions here. 
# Supposedly it Inhibits osteoclast formation.

# Because all of the top mean variance have very small expressions. I'm 
# going to only look at mean variances where the expressions are above the 
# 25th percentile of avg expression


#---------------------------------------------------------------------------
# Investigating mean - variance for genes above 25% percentile for avg expression
#---------------------------------------------------------------------------
first_quartile <- summary(row_mean)[2]


summary_matrix_high <- summary_matrix[row_mean >=  first_quartile,]
summary_matrix_high

df_summary_matrix_high <- data.frame(summary_matrix_high)
ggplot(df_summary_matrix_high) + 
  geom_histogram(mapping = aes(row_rel_sd))

summary_matrix_high <- summary_matrix_high [names(sort(summary_matrix_high[,4], decreasing = TRUE)), ]

top_5_percent <- summary_matrix_high[ 1 : (0.05*nrow(summary_matrix_high)),]
rownames(top_5_percent) <- rownames(summary_matrix_high[1 : (0.05*nrow(summary_matrix_high)),])

# top_5_percent has the top 5% highest mean_variances

#---------------------------------------------------------------------------
# Investigating high mean variances
#---------------------------------------------------------------------------

#Bpifa1 : gram-negative defence in airways
# I'd expect it to be present later on when nasal passages have formed

tfs <- rownames(top_5_percent[1:50,])
tfs
top_5_percent[50:60,]

plot_list <- lapply(tfs, function(x) {
  
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x)
})
names(plot_list) <- tfs

# test_intervals <- list(c(1:10),
#                     c(11:20),
#                     c(21:30),
#                     c(31:40),
#                     c(41:50),
#                     c(51:60),
#                     c(61:70),
#                     c(71:80),
#                     c(81:90),
#                     c(91:100))
# n <- 1
# while (n+10 <100) {
#   y <- n+10
# 
#   tfs <- rownames(top_5_percent[n:y,])
#   
#   plot_list <- lapply(tfs, function(x) {
#     
#     df <- data.frame(Counts = count_matrix[x, ],
#                      Dev_stage = meta_data$dev_stage)
#     
#     df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
#     
#     ggplot(df, aes(y = Counts, x = Dev_stage)) +
#       geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
#       ggtitle(x)
#   })
#   names(plot_list) <- tfs
#   
#   view(plot_grid(plotlist = plot_list, nrow = 2, ncol = 5))
#   
#   n <- n+10
# }
# I'm not sure why the plots arn't showing up in this while loop

 #It seems like most of these just have high variances because they have outliers
# Likely need to find a way to filter out genes that don't have outliers. But 
# you can't really determine which are outliers because there's only n=2 for each 
# dev stage. 


# VIP looks like it makes sense as its a gstrointestinal gene. So it makes sense
# that it is expressed so highly in P0

# I'm going to 


#---------------------------------------------------------------------------
# Exploring cell types and organism parts
#--------------------------------------------------------------------------

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5754028/
# hippocampus dev begins at E11

# http://hippocampome.org/php/markers.php
# Compile positive markers

tfs <- c("Calb1","Nrgn","Actn2")
plot_list <- lapply(tfs, function(x) {
  
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)

# For Calb1 and Nrgn, you can actually see that they increse at 11.5
# Calb1 increase significantly at E11.5, Whereas Nrgn starts increasing at 11.5
# So you can sortof see CA1 Pyramidal formation

# Compile negative markers
tfs <- c("Calb2","Pvalb","Cck","Sst","Htr3a", "Reln")
plot_list <- lapply(tfs, function(x) {
  
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)

count_matrix["Calb2",]

# Unfortunately all of these negative markers don't seem to be decreasing over time.
# its possible these markers may be negative for hippocampus, but 
# they may be positive for other things.

# I realize that i'm stupid for choosing hippocakmpus because 
# our samples are all forebrain. oops


#---------------------------------------------------------------------------
# Exploring cell types of the forebrain
#--------------------------------------------------------------------------

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5095704/

#Basal forbrain has 3 common cell types

#cholinergic
#glutamatergic
#GABAergic


# --- Exploring Cholinergic Neuons

#https://www.abcam.com/neuroscience/cholinergic-neuron-markers-and-their-functions

#Markers: 
#Choline acetyltransferase (ChAT)
#Vesicular acetylcholine transporter(VACht) 
# Acetylcholinesterase

tfs <- c("Chat","Slc18a3","Ache")
plot_list <- lapply(tfs, function(x) {
  
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)

# Nice trend increase in cholinergic neuron markers




# --- Exploring glutamatergic Neurons
# Glutamate is the neurotransmitter

#Markers: https://www.abcam.com/neuroscience/glutamatergic-neuron-markers-and-their-functions

#vGluT2: "Slc17a6"
#vG1uT1: Slc17a7
#NMDAR1: Grin1
#NMDAR2B: Grin2b

tfs <- c("Slc17a6", "Slc17a7", "Grin2b", "Grin1")

plot_list <- lapply(tfs, function(x) {
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)


# --- Exploring GABAergic Neurons


#Markers:https://www.abcam.com/neuroscience/gabaergic-neuron-markers-and-their-functions

#
#Gat1: Slc6a1
#Gaba receptor 1: Gabbr1
#Gaba receptor 2: Gabbr2
# Glutamate decarboxylase isoforms65 : Gad2
#Glutamate decarboxylase isoforms 67: Gad1

tfs <- c("Slc6a1", "Gabbr1", "Gabbr2", "Gad2", "Gad1")

plot_list <- lapply(tfs, function(x) {
  df <- data.frame(Counts = count_matrix[x, ],
                   Dev_stage = meta_data$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)
