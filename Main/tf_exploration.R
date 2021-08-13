
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
# Log 10+1 normalizing 
#---------------------------------------------------------------------------
count_matrix <- log10(count_matrix+1)
avg_count_matrix <- log10(avg_count_matrix+1)


#---------------------------------------------------------------------------
# Create place for new matrix exploration graphs 
#---------------------------------------------------------------------------
dir.create(path = "Data/TF_Exploration")
dir.create(path = "Data/TF_Exploration/General_Graphs")

#---------------------------------------------------------------------------
# Exploring TFs requried for brain morphology
#--------------------------------------------------------------------------

v_tfs_req <- c("Chrd", "Nog", "Hesx1")

df_tfs_req <- as.data.frame(avg_count_matrix[v_tfs_req,])%>%
  mutate("gene" = rownames(avg_count_matrix[v_tfs_req,]))

df_tfs_req <- melt(df_tfs_req, id.vars <- c("gene"))

#want ggplot x = dev stage, y = expression, facet wrap on gene


pl_tfs_req <- ggplot(df_tfs_req) +
  geom_point(mapping= aes(x = variable,
                          y = value), fill = "steelblue")+
  facet_wrap( ~gene) + 
  xlab("Developmental Stage") +
  ylab ("Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  


pl_tfs_req


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

tfs_chol <- c("Chat","Slc18a3","Ache")

df_tfs_chol <- as.data.frame(avg_count_matrix[tfs_chol,])%>%
  mutate("gene" = rownames(avg_count_matrix[tfs_chol,]))

df_tfs_chol <- melt(df_tfs_chol, id.vars <- c("gene"))

pl_tfs_chol<- ggplot(df_tfs_chol) +
  geom_point(mapping= aes(x = variable,
                          y = value), fill = "steelblue")+
  facet_wrap( ~gene) + 
  xlab("Developmental Stage") +
  ylab ("Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pl_tfs_chol

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

tfs_gaba <- c("Slc6a1", "Gabbr1", "Gabbr2", "Gad2", "Gad1")

df_tfs_gaba <- as.data.frame(avg_count_matrix[tfs_gaba,])%>%
  mutate("gene" = rownames(avg_count_matrix[tfs_gaba,]))

df_tfs_gaba <- melt(df_tfs_gaba, id.vars <- c("gene"))

pl_tfs_gaba<- ggplot(df_tfs_gaba) +
  geom_point(mapping= aes(x = variable,
                          y = value), fill = "steelblue")+
  facet_wrap( ~gene) + 
  xlab("Developmental Stage") +
  ylab ("Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pl_tfs_gaba




#---------------------------------------------------------------------------
# Subsetting Count matrix for only Transcription Factors
#---------------------------------------------------------------------------

all_tfs <- read.delim(file = "mouse_tfs.tsv", 
                  sep = "\t",
                  stringsAsFactors = FALSE)

tf_count_matrix <- count_matrix[rownames(count_matrix)%in% all_tfs$Symbol ,]

df_tf_expressed <- data.frame(numbers = c(
                                          nrow(count_matrix),
                                          nrow(all_tfs),
                                          nrow(tf_count_matrix)),
  
                              names = factor(c("PC Genes",
                                        "All Tfs",
                                        "Expressed Tfs"),
                                        levels = c("PC Genes",
                                                    "All Tfs",
                                                    "Expressed Tfs")))

pl_tf_expressed <- ggplot(df_tf_expressed) + 
  geom_col(mapping = aes( x = names, y = numbers), fill = "steelblue") +
  ggtitle("Subsetting Transcription Factors") +
  ylab("Count")

pl_tf_expressed


dir.create("Data/TF_Exploration/Alex_Tfs")
ggsave(filename = "Data/TF_Exploration/Alex_Tfs/Tf_mapping.png")

#---------------------------------------------------------------------------
# Summary stats for Transcription Factor Subset
#---------------------------------------------------------------------------
tf_row_sum <- rowSums(tf_count_matrix)
tf_row_mean <- rowMeans(tf_count_matrix)
tf_row_stdev <- rowSD(tf_count_matrix)
tf_row_rel_sd <- tf_row_stdev/tf_row_mean

tf_summary_matrix <- matrix(c(tf_row_sum,
                           tf_row_mean,
                           tf_row_stdev,
                           tf_row_rel_sd),
                         nrow = length(tf_row_sum),
                         ncol = 4)
colnames(tf_summary_matrix) <-  c("tf_row_sum",
                               "tf_row_mean",
                               "tf_row_stdev",
                               "tf_row_rel_sd")
rownames(tf_summary_matrix) <-  names(tf_row_sum)

tf_summary_matrix

summary_tf <- summary(tf_summary_matrix[,"tf_row_mean"])

pl_tf_avg <- ggplot(as.data.frame(tf_summary_matrix)) +
  geom_histogram(mapping = aes(tf_row_mean), fill = "steelblue")+
  xlab("Expression (log10+1 normalized)") +
  ylab("Count")+
  ggtitle("Expression of Transcription Factors") + 
  geom_vline(xintercept = summary_tf["Mean"], color = "red")
  
ggsave(filename = "Data/TF_Exploration/Alex_Tfs/avg_tfs.png")
  
#---------------------------------------------------------------------------
# Exploring Transcription Factors
#---------------------------------------------------------------------------

view(tf_summary_matrix[ order(tf_summary_matrix[ , "tf_row_mean"], 
                              na.last = TRUE,
                              decreasing = TRUE) ,  ])
#Highest mean is "Hmgb1"



#---------------------------------------------------------------------------
# Exploring Alex's 8 Tfs
#---------------------------------------------------------------------------

tfs_alex <- c("Ascl1",
              "Tcf4", 
              "Neurod1",
              "Runx1",
              "Mecp2", 
              "Mef2c",
              "Hes1", 
              "Pax6")

tf_matr_alex <- tf_count_matrix[tfs_alex,]

tf_matr_avg_alex <- avg_count_matrix[tfs_alex,]

df_tfs_alex <- as.data.frame(tf_matr_avg_alex)%>%
  mutate("gene" =rownames(tf_matr_avg_alex))

df_tfs_alex <- melt(df_tfs_alex, id.vars <- c("gene"))

pl_tfs_alex<- ggplot(df_tfs_alex) +
  geom_point(mapping= aes(x = variable,
                          y = value), fill = "steelblue")+
  facet_wrap( ~gene) + 
  xlab("Developmental Stage") +
  ylab ("Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pl_tfs_alex

ggsave(filename = "Data/TF_Exploration/Alex_Tfs/Alex8.png")




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



