library(tidyverse)
library(preprocessCore)
library(cowplot)

meta <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE)
count_mat <- readRDS(file = "Data/pc_count_matrix.rds")
avg_count_mat <- readRDS(file = "Data/avg_pc_count_matrix.rds")
diff_count_mat <- readRDS(file = "Data/diff_pc_count_matrix.rds")


# Process: quantile norm (log10+1) counts

count_mat <- log10(count_mat+1)
proc_mat <- normalize.quantiles(count_mat)


# look at the distribution of gene counts before and after norm
# data is already very tight - quantile norm has minimal effect

plot(density(count_mat)) # first sample

for (i in 2:ncol(count_mat)) {  # rest of the samples
  lines(density(count_mat[, i]), col = "black")
}

for (i in 1:ncol(proc_mat)) {  # then all of the normalized samples
  lines(density(proc_mat[, i]), col = "red")
}

# Sum of counts across samples - use to look for no/low expressed
gene_counts <- rowSums(count_mat)
summary(gene_counts)
plot(density(gene_counts))
hist(gene_counts, breaks = 100)
sum(gene_counts == 0)

# look at the gene averages across all samples

gene_avg <- rowMeans(count_mat)
summary(gene_avg)
head(sort(gene_avg, decreasing = TRUE), 20)

# look at the standard deviation across all samples

gene_sd <- apply(count_mat, 1, sd)
summary(gene_sd)
head(sort(gene_sd, decreasing = TRUE), 20)
plot(count_mat["Neurod6",], type = "line")

# mean variance relationship

plot(gene_avg, gene_sd)


# Looking at TFs of interest

tfs <- c("Ascl1", "Tcf4", "Neurod1", "Runx1", "Mecp2", "Mef2c", "Hes1", "Pax6")

plot_list <- lapply(tfs, function(x) {
  
  df <- data.frame(Counts = count_mat[x, ],
                   Dev_stage = meta$dev_stage)
  
  df$Dev_stage <- factor(df$Dev_stage, levels = unique(df$Dev_stage))
  
  ggplot(df, aes(y = Counts, x = Dev_stage)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "royalblue") +
    ggtitle(x) +
    theme_minimal() +
    theme(axis.title.x = element_blank())
})

names(plot_list) <- tfs

plot_list$Ascl1
plot_grid(plotlist = plot_list, nrow = 2, ncol = 4)
