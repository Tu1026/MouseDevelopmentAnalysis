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
library(edgeR)


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
count_matrix <- log10(count_matrix+1)
avg_count_matrix <- log10(avg_count_matrix+1)


#---------------------------------------------------------------------------
# Create place for new matrix exploration graphs 
#---------------------------------------------------------------------------
dir.create(path = "Data/Matrix_Exploration")
dir.create(path = "Data/Matrix_Exploration/General_Graphs")

#---------------------------------------------------------------------------
# Identifying differentially expressed genes between replicates
#---------------------------------------------------------------------------
#Pass
#---------------------------------------------------------------------------
# Identifying most expressed genes
#---------------------------------------------------------------------------

# Here we will identify the genes that are most expressed in each sample

# --- Looking at a single sample first. Lets look at E14.5-2

# Graph histogram For E14.5-2

E14 <- as.data.frame(count_matrix[,14])

general_hist <- ggplot(data = E14) +
  geom_histogram(mapping = aes(count_matrix[,14]), fill = "steelblue")+
  xlab("Expression (log10+1 normalized)") + 
  ylab("Count") + 
  ggtitle("Expression Histogram of E14.5-2")

ggsave(filename = "Data/Matrix_Exploration/General_Graphs/histogram.png",
       general_hist)


E14_most_expr <- arrange(E14, desc(E14))
colnames(E14_most_expr) <- c("Log10+1 Expression")

head(E14_most_expr,5)

# --- Doing the same for all dev stages
l_most_expr_genes <- list()

get_most_expr<- function(column) {
  my_df <- data.frame("count" = column)
  arranged_df <- arrange(my_df, desc(my_df))
  return(top_n(arranged_df, n = 10))
}

most_expr_genes <- apply(count_matrix, MARGIN = 2, FUN = get_most_expr)
head(most_expr_genes, 5)

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
head(never_expressed_genes, 5) 

#---------------------------------------------------------------------------
# Summarizing the count_matrix
#---------------------------------------------------------------------------
#Get summary stats for the count_matrix and make a boxplot

summary <- summary(count_matrix[,])

count_boxplot <- ggplot(stack(as.data.frame(count_matrix)))+
  geom_boxplot(mapping = aes(x = ind, y = values))+
  xlab("Sample") + 
  ylab("Expression (log10+1 normalized)") + 
  ggtitle("Boxplot of Samples") +
  scale_x_discrete(labels = c("E10.5-1",
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
  
count_boxplot

ggsave(filename = "Data/Matrix_Exploration/General_Graphs/count_boxplot.png",
       count_boxplot)





#---------------------------------------------------------------------------
# Genes only expressed at P0
#---------------------------------------------------------------------------
#pull an example just because

never_expr_p0 <- count_matrix[  count_matrix[ , 1] <= log10(1) &
                           count_matrix[ , 2] <= log10(1) & 
                           count_matrix[ , 3] <= log10(1) &
                           count_matrix[ , 4] <= log10(1) &
                           count_matrix[ , 5] <= log10(1) &
                           count_matrix[ , 6] <= log10(1) &
                           count_matrix[ , 7] <= log10(1) &
                           count_matrix[ , 8] <= log10(1) &
                          count_matrix[ , 9] <= log10(1) & 
                          count_matrix[ , 10] <= log10(1) &
                          count_matrix[ , 11] <= log10(1) &
                          count_matrix[ , 12] <= log10(1) &
                          count_matrix[ , 13] <= log10(1) &
                          count_matrix[ , 14] <= log10(1) , ]

#Retnla is the best case

count_matrix["Retnla",] # it is clearly only expressed in P0

retnla<- ggplot(as.data.frame(avg_count_matrix["Retnla",])) +
  geom_point(mapping= aes(x = colnames(avg_count_matrix),
                          y = avg_count_matrix["Retnla",]),
             fill = "steelblue")+
  xlab("Developmental Stage") +
  ylab ("Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Expression of Retn1a")


retnla
ggsave(filename = "Data/Matrix_Exploration/General_Graphs/Retnla.png")

#---------------------------------------------------------------------------
# Summarizing Gene Expression across rows and constructing matrix
#---------------------------------------------------------------------------

# For convenience, I'm going to make a little matrix that has all of the row summary
# information in it .

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
# sTdev convering at 4
#---------------------------------------------------------------------------

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
# Genes that are expressed
#---------------------------------------------------------------------------

expressed_genes <- summary_matrix [row_sum != 0,]
expressed_genes

stopifnot(nrow(never_expressed_genes) + nrow(expressed_genes) == nrow(summary_matrix))

# Genes that are ordered by row_sum
expressed_genes[order(expressed_genes[,"row_sum"]), ]

# Proof that it works. Values are the same
expressed_genes["Dppa2",]

summary(expressed_genes[,"row_sum"])

nrow(expressed_genes)


df_expressed_genes <- data.frame("Sum_Expression" = expressed_genes[,"row_sum"],
                                 "Genes" = rownames(expressed_genes),
                                 "RSD_Expression" = expressed_genes[,"row_rel_sd"],
                                 "SD_Expression" = expressed_genes[, "row_stdev"])
ggplot(df_expressed_genes) + 
  geom_histogram(mapping = aes(Sum_Expression), fill = "red") +
  geom_vline(xintercept = c(summary(expressed_genes[,"row_sum"]))) +
  geom_text(aes(x=summary(expressed_genes[,"row_sum"])[3],
                label="Median",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11)) +
  geom_text(aes(x=summary(expressed_genes[,"row_sum"])[4],
                label="Mean",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11))+
  labs(title = "Distribution of row_sum for expressed genes")

ggplot(df_expressed_genes) + 
  geom_histogram(mapping = aes(RSD_Expression), fill = "red", color = "black") +
  geom_vline(xintercept = c(summary(expressed_genes[,"row_rel_sd"]))) +
  geom_text(aes(x=summary(expressed_genes[,"row_rel_sd"])[3],
                label="Median",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11)) +
  geom_text(aes(x=summary(expressed_genes[,"row_rel_sd"])[4],
                label="Mean",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11))+
  labs(title = "Distribution of row_rel_sd for expressed genes")

#---------------------------------------------------------------------------
# Why values stack up at 4 for relative standard deviation
#---------------------------------------------------------------------------

# Case example is "Cpn2" or "Cdcp3"

count_matrix["Cdcp3",]
count_matrix["Cpn2",]

origonal_count <- (10 ** count_matrix["Cdcp3",]["ENCFF918QNL"]) - 1

# I'm presuming that 0.01 was the lowest possible value for their TPM

# so the values of 4 are when one instance of the loest possible epxression happens


count_matrix["Pigr",]

# We can see that rel st values of 2.73 correspond to when there are 2 instances 
# of the loest possible expression


#---------------------------------------------------------------------------
# mean variance plot for expressed genes
#---------------------------------------------------------------------------

plot(expressed_genes[,"row_mean"], expressed_genes[,"row_stdev"])


#---------------------------------------------------------------------------
# Removing lowly expressed Genes with high variance
#---------------------------------------------------------------------------

# I'm going to say that in order to be expressed, the row_sum of a gene must be
# greater than 1.000. Using the log transformed values this is around 0.30103
expressed_genes[ order (expressed_genes[,"row_sum"]) , ]

high_expressed_genes <- expressed_genes[ expressed_genes[ ,"row_sum"] > log10(2), ]
high_expressed_genes[ order (high_expressed_genes[,"row_sum"]) , ]

df_high_expressed_genes <- data.frame(
  "Sum_Expression" = high_expressed_genes[,"row_sum"],
  "RSD_Expression" = high_expressed_genes[,"row_rel_sd"])


ggplot(df_high_expressed_genes) + 
  geom_histogram(mapping = aes(Sum_Expression), fill = "red") +
  geom_vline(xintercept = c(summary(high_expressed_genes[,"row_sum"]))) +
  geom_text(aes(x=summary(high_expressed_genes[,"row_sum"])[3],
                label="Median",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11)) +
  geom_text(aes(x=summary(high_expressed_genes[,"row_sum"])[4],
                label="Mean",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11))+
  labs(title = "Distribution of row_sum for highly_ expressed genes")

ggplot(df_high_expressed_genes) + 
  geom_histogram(mapping = aes(RSD_Expression), fill = "red") +
  geom_vline(xintercept = c(summary(high_expressed_genes[,"row_rel_sd"]))) +
  geom_text(aes(x=summary(high_expressed_genes[,"row_rel_sd"])[3],
                label="Median",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11)) +
  geom_text(aes(x=summary(high_expressed_genes[,"row_rel_sd"])[4],
                label="Mean",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11))+
  labs(title = "Distribution of row_rel_sd for highly_expressed genes")


# Its important to note that the distribution for the highly expressed genes
# doesn't have that stack up at 4 For RSD. 

# Notably, all of the higehst mean variance genes have really low expression


#---------------------------------------------------------------------------
# Investigatging Some random genes
#---------------------------------------------------------------------------
# I've chosen this gene artibrarily. Just to see what it looks like over time

# Some lowly expressed genes

expressed_genes[ order (expressed_genes[,"row_sum"])  , ]

plot(count_matrix["Olfr183",])
plot(count_matrix["Serpinb3b",])

# Some lowly expressed Genes from the bottom of the high_expressed_genes
high_expressed_genes[ order (high_expressed_genes[,"row_sum"]) , ]
summary(high_expressed_genes)

plot(count_matrix["Trim43a",])
plot(count_matrix["Btla",])

# Some highly expressed Genes from between 1rst and 2nd percentile
high_expressed_genes[ high_expressed_genes[,"row_sum"] > 5.27, ]
plot(count_matrix["Itga3",])

# Some highly expressed Genes from between 2nd and 3rd percentile
plot(count_matrix["Kdelr1",])

# Some highly expressed Genes from between 3rd and 4th percentile
high_expressed_genes[ high_expressed_genes[,"row_sum"] > 6.67, ]
plot(count_matrix["Ap1s1",])

# Most Highly Expressed Genes
high_expressed_genes[ order (high_expressed_genes[,"row_sum"], decreasing = TRUE ) , ]

plot(count_matrix["Tuba1a",])
plot(count_matrix["Stmn1",])
plot(count_matrix["Ptges3",])


#---------------------------------------------------------------------------
# Investigating mean - variance for genes above 25% percentile for avg expression
#---------------------------------------------------------------------------
first_quartile <- summary(row_mean)[2]

# Double check all of the first quartile's and see if any are small expression to high expression
# Plot the before

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

# TODO Make slide about how yea these are increasing and thats indicative of growth/increasing abundance of cell that expresses
# TODO Be able to talk about the reliability of this data


#---------------------------------------------------------------------------
# Exploring Tf factors
#---------------------------------------------------------------------------

# explore their expression changes

tfs <- read.delim(file = "mouse_tfs.tsv", 
                  sep = "\t")

tf_count_matrix <- count_matrix[tfs$Symbol,]

tf_summary_matrix <- summary_matrix[tfs$Symbol,]

stopifnot(nrow(tfs) == nrow(tf_count_matrix))
stopifnot(nrow(tfs) ==  nrow(tf_summary_matrix))

# --- Exploring expression change

df_tf_count <- data.frame("TF_Sum" = tf_summary_matrix[,"row_sum"])

ggplot(df_tf_count) + 
  geom_histogram(mapping = aes(TF_Sum), fill = "red") +
  geom_vline(xintercept = c(summary(tf_summary_matrix[,"row_sum"]))) +
  geom_text(aes(x=summary(tf_summary_matrix[,"row_sum"])[3],
                label="Median",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11)) +
  geom_text(aes(x=summary(tf_summary_matrix[,"row_sum"])[4],
                label="Mean",
                y=20),
            colour="blue", angle=90, vjust = 1.2, text=element_text(size=11))+
  labs(title = "Distribution of row_sum for tfs expressed genes")

# Its a very similar expression pattern to the overall genes. Not surprisingly.


# --- Correlating TF factors to eachother

# Generate a corr matrix for each TF across the time samples

rot_count_matrix <- t(count_matrix)

# TODO FORALEX: Is this what you want? It takes a long time. How do I pick the strongest correlations?

rot_pear_rcorr_object <- rcorr(rot_count_matrix[,1:100], type = "pearson")
rot_pear_matrix <- rot_pear_rcorr_object$r
rot_pear_matrix_pvalues <- rot_pear_rcorr_object$P

rot_spear_heatmap <- heatmap(x = rot_pear_matrix,
                             sym = TRUE,
                             Rowv = NA,
                             Colv = NA, 
                             revC= TRUE)

# ###HELP###
# for (colname in colnames(rot_pear_matrix)) {
#   # browser()
#   
#   bool <- rot_pear_matrix[,colname] > 0.9 & rot_pear_matrix[,colname] != 1.00
#   
#   rot_pear_matrix[,bool]
# }
# view(rot_pear_matrix)

# ------ Correlate TF factors to other strongly positive correlated genes

# Lets look at Gnai3 and Cdc45 because they are supposed to have 0.9 correlation

plot(count_matrix["Gnai3", ])
plot(count_matrix["Cdc45", ])

#You can see that they are very closely positive correlated

plot(count_matrix["Gnai3", ],
     count_matrix["Cdc45", ])

# Pretty correlated


# ------ Correlate TF factors to other strongly negative correlated genes

# Lets look at Gnai3 and Clec10a because they are supposed to have -0.9 correlation

plot(count_matrix["Gnai3", ])
plot(count_matrix["Clec10a", ])

#You can see that they are very strongly negatively correlated

plot(count_matrix["Gnai3", ],
     count_matrix["Clec10a", ])

# Pretty correlated


# ------ Correlate TF factors to non-correlated genes

# Lets look at Gnai3 and Gpa33 because they are supposed to have ~0 correlation
# We are asuming GPa33 is normally distributed, which it prob isnt because its variation is so low

plot(count_matrix["Gnai3", ])
plot(count_matrix["Gpa33", ])

#You can see that they are very weakly correlated

plot(count_matrix["Gnai3", ],
     count_matrix["Gpa33", ])

# Pretty non- correlated


# TODO
# Generate matrices forother parts

#TODO design matrix

# TODO DEA

# TODO make plots nier and begin making good slides

#### 
#Remove low counts and save that as count_matrix_pro