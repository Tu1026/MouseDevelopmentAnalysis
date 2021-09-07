# install.packages("infer")
# 
library(infer)
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

count_matrix <- readRDS(file = "Data/pc_count_matrix0slog2.rds")
avg_count_matrix <- readRDS(file = "Data/avg_pc_count_matrix0slog2.rds")
meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")

#---------------------------------------------------------------------------
# Removing 0s
#---------------------------------------------------------------------------
# count_matrix <- count_matrix[ count_matrix[ , 1] > 0 |
#                               count_matrix[ , 2] > 0 | 
#                               count_matrix[ , 3] > 0 |
#                               count_matrix[ , 4] > 0 |
#                               count_matrix[ , 5] > 0 |
#                               count_matrix[ , 6] > 0 |
#                               count_matrix[ , 7] > 0 |
#                               count_matrix[ , 8] > 0 |
#                               count_matrix[ , 9] > 0 | 
#                               count_matrix[ , 10] > 0 |
#                               count_matrix[ , 11] > 0 |
#                               count_matrix[ , 12] > 0 |
#                               count_matrix[ , 13] > 0 |
#                               count_matrix[ , 14] > 0 |
#                               count_matrix[ , 15] > 0 |
#                               count_matrix[ , 16] > 0 , ]
# 
# avg_count_matrix <- avg_count_matrix[ avg_count_matrix[ , 1] > 0 |
#                                 avg_count_matrix[ , 2] > 0 | 
#                                 avg_count_matrix[ , 3] > 0 |
#                                 avg_count_matrix[ , 4] > 0 |
#                                 avg_count_matrix[ , 5] > 0 |
#                                 avg_count_matrix[ , 6] > 0 |
#                                 avg_count_matrix[ , 7] > 0 |
#                                 avg_count_matrix[ , 8] > 0 , ]

#---------------------------------------------------------------------------
# Log2 normalizing 
#---------------------------------------------------------------------------
# count_matrix <- log2(count_matrix+1)
# avg_count_matrix <- log2(avg_count_matrix+1)


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
  ylab ("Expression (Log2 Transformed") +
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
  ylab ("Expression (Log2 Transformed)") +
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
  ylab ("Expression (Log2 Transformed)") +
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
                                          nrow(tf_count_matrix)),
                              names = factor(c("PC Genes",
                                            "Expressed Tfs"),
                                            levels = c("PC Genes",
                                                       "Expressed Tfs")
                                            ))



pl_tf_expressed <- ggplot(df_tf_expressed) + 
  geom_col(mapping = aes( x = names, y = numbers), fill = "steelblue") +
  ggtitle("Subset of Transcription Factors") +
  ylab("Count")+
  xlab("")

pl_tf_expressed


dir.create("Data/TF_Exploration/Alex_Tfs", showWarnings = FALSE)
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
  xlab("Mean Expression (Log2 Transformed)") +
  ylab("Count")+
  ggtitle("Mean Expression of Transcription Factors") + 
  geom_vline(xintercept = summary_tf["Mean"], color = "red")

tf_summary_matrix[]
 
pl_tf_avg 
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
  facet_wrap( ~gene,
              scales = "free") + 
  xlab("Developmental Stage") +
  ylab ("Expression (Log2 Transformed)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pl_tfs_alex



ggsave(filename = "Data/TF_Exploration/Alex_Tfs/Alex8.png")

#--Get Ruxn1 graph
runx1 <- df_tfs_alex %>%
  filter(gene == "Runx1")

pl_runx1<- ggplot(runx1) +
  geom_point(mapping= aes(x = variable,
                          y = value), fill = "steelblue")+
  xlab("Developmental Stage") +
  ylab ("Expression (Log2 Transformed)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(title = "Average Runx1 Expression")


# --- Exploring expression change


ggplot(df_tf_count) + 
  geom_histogram(mapping = aes(TF_Sum), fill = "red") +
  geom_vline(xintercept = c(summary(tf_summary_matrix[,"tf_row_sum"]))) +
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


#---------------------------------------------------------------------------
# The Reproducibility of Mecp2
#---------------------------------------------------------------------------
Mecp2 <- data.frame(Expression = tf_count_matrix["Mecp2",] ,
                    dev_stage = factor(meta_data$dev_stage),
                    group = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8))
ggplot(Mecp2) +
  geom_point(mapping = aes(x = dev_stage, y = Expression))


Mecp2

#---------------------------------------------------------------------------
# Isolating Runx1 graph for powerpoint
#---------------------------------------------------------------------------
Runx1_avg <- data.frame(avg = avg_count_matrix["Runx1",] ,
                        dev_stage = factor(unique(meta_data$dev_stage)))
Runx1_avg$avg
Runx1_avg$dev_stage


Runx1_graph_avg<- ggplot(Runx1_avg) +
  geom_point(mapping= aes(x = dev_stage,
                          y = avg),
             fill = "steelblue") +
  xlab("Developmental Stage") +
  ylab ("Avg Expression (log10+1 normalized)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs( title = "Avg Runx1 Expression")

Runx1_graph_avg
ggsave(filename = "Data/TF_Exploration/Alex_Tfs/Runx1.png")



                                        
#---------------------------------------------------------------------------
# Finding Genes highly correlated with Ruxn1 expression
#---------------------------------------------------------------------------
pear_rcorr_object <- rcorr(t(tf_count_matrix), type = "pearson")
pear_matrix <- pear_rcorr_object$r
pear_matrix_pvalues <- pear_rcorr_object$P

pear_matrix[1:10,]

# pear_matrix[(upper.tri(pear_matrix))] <- NA # Remove the entire upper triangle

pear_heatmap <- pheatmap(pear_matrix, 
                         main = "Runx1 Correlations",
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         breaks = rev(c(1.0,0.99999999, 0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6)),
                         color = rev(brewer.pal(n = 10, name ="RdYlBu")),
                         legend_breaks = c(1.0,0.9,0.8,0.7,0.6),
                         na_col = "blue",
                         border_color = "white",
                         show_rownames = FALSE,
                         show_colnames = FALSE
)
pear_heatmap

Runx1_pears <- data.frame(Pearson = pear_matrix["Runx1",])

Runx1_pears <- Runx1_pears %>%
  arrange(Pearson)
Runx1_pears <- Runx1_pears %>%
  mutate(gene = rownames(Runx1_pears))

#---------------------------------------------------------------------------
# Highest and lowest pearson correlations to RUnx1
#---------------------------------------------------------------------------
Runx1_tail <- tail(Runx1_pears, 10)
Runx1_head <- head(Runx1_pears, 10)


#---------------------------------------------------------------------------
# Chipseq and Pearson Bringint the two together
#---------------------------------------------------------------------------

chip <- readRDS(file = "Data/mouse_mean_bind.RDS")

Runx1_chip <- data.frame(Binding = chip$Runx1)
any(duplicated(rownames(Runx1_chip)))

Runx1_chip_keep <- Runx1_chip %>%
  filter(rownames(Runx1_chip) %in% rownames(Runx1_pears))
Runx1_chip_keep <-  Runx1_chip_keep %>%
  mutate(gene = rownames(Runx1_chip_keep))

Runx1_pears_chip <- inner_join ( Runx1_pears , Runx1_chip_keep, by = "gene")

df <- data.frame(Runx1_pears_chip)

ggplot(df) +
  geom_point(mapping = aes(x = Binding, y = Pearson), colour = "steelblue") +
  labs( title = "Pearson Correlation vs Mean Binding Score", y= "Pearson Correlation",x = "Mean Binding Score")

ggplot(df) 

#---------------------------------------------------------------------------
# Adding on the n_times ??
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# Investigating High correalted Runx1 genes and their chipseq scores
#---------------------------------------------------------------------------

Runx1_upreg <- filter(Runx1_pears_chip, gene %in% Runx1_tail$gene) %>%
  arrange(by_group = desc(Pearson))
Runx1_downreg <- filter(Runx1_pears_chip, gene %in% Runx1_head$gene) %>%
  arrange(by_group = desc(Pearson))



















#---------------------------------------------------------------------------
# Try to get P-value for our chipseq scores
#---------------------------------------------------------------------------

# Remove low Runx1chip data
Runx1_pro <- filter(Runx1_pears_chip, !(Binding <0.1))
nrow(Runx1_pro)
nrow(Runx1_pears_chip)

hist(Runx1_pro$Binding)
skim(Runx1_pears_chip$Binding)
hist(Runx1_pears_chip$Binding)


ggplot(Runx1_pro) + 
  geom_histogram(mapping = aes(Binding), fill = "steelblue", color = "black")+
  labs(x = "Mean Binding Score", y = "Count", title = "Chip-seq Binding Scores for Runx1 Pulldown")
# don-t need to remove low data, fairly normally distributed

#---------------------------------------------------------------------------
# Try to get P-value for our chipseq scores
#--------------------------------------------------------------------------

# Do z-tests and merge
qqnorm(Runx1_upreg$Binding)
qqline(Runx1_upreg$Binding)
# not normally distributed enough?


# T-test (even though its wrong, we should be using z-test)
p_value <- c()
for (n_row in 1:nrow(Runx1_upreg)) { 
  my_t <- t.test(x = Runx1_pears_chip$Binding, mu = Runx1_upreg[n_row,"Binding"])
  p_value <- append(p_value, my_t$p.value)
}
p_value

Runx1_upreg <- cbind(Runx1_upreg,p_value)


#Infer package hypothesis testing
mean(Runx1_pro$Binding)

upreg_null <- Runx1_pro %>%
  specify(response = Binding)%>%
  hypothesise(null = "point", mu = mean(Runx1_pro$Binding)) %>%
  generate(reps = 10000, type = "bootstrap") %>%
  calculate(stat = "mean")
visualize(upreg_null)+
  shade_p_value(Runx1_upreg[3,"Binding"], direction = "left")
#oh god

