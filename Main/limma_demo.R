# Demonstration of using limma to find differentially expressed genes

library(limma)
library(tidyverse)
library(edgeR)
library(tximport)

# load data
meta <- read.delim("Main/Data/complete_meta_data.tsv", stringsAsFactors = FALSE)
pc <- read.delim("Main/Data/pc_sub.tsv", stringsAsFactors = FALSE)
files <- paste0("Main/Data/Count_tables/", meta$id, ".tsv")
counts <- tximport(files, type = "rsem", txIn = FALSE)
count_mat <- counts$counts

# name matrix and only keep protein coding genes
rownames(count_mat) <- str_replace(rownames(count_mat), "\\.[:digit:]*$", "")
colnames(count_mat) <- meta$id
count_mat <- count_mat[rownames(count_mat) %in% pc$Gene_ID, ]

# order by tissue and dev stage
meta <- meta %>% group_by(tissue_type) %>% arrange(dev_stage, .by_group = TRUE)
count_mat <- count_mat[, meta$id]
stopifnot(identical(colnames(count_mat), meta$id))

# build design matrix
design <- model.matrix(~ 0 + meta$tissue_type + meta$dev_stage)

# edgeR: remove lowly expressed, scale normalize, and the log2 CPM transform
dge <- DGEList(count_mat)
keep <- filterByExpr(dge, design)
dge <- dge[keep, ]

# gene counts of removed
removed <- rowSums(count_mat[!keep,])
summary(removed)
count_mat[which.max(removed), ]

# scale normalization using TMM method (default)
dge <- calcNormFactors(dge)

# convert counts to log2CPM 
cpm_counts <- cpm(dge, log = TRUE, prior.count = 3)

# fit model with limma eBayes trend 
fit <- lmFit(cpm_counts, design)
fit <- eBayes(fit, trend=TRUE)

# pulling specific contrasts
test_age <- topTable(fit, number = Inf, coef = "meta$dev_stageP0")
test_tissue <- topTable(fit, number = Inf, coef = "meta$tissue_typethymus")

plot(cpm_counts[rownames(test_age)[1], ], col = ifelse(meta$dev_stage == "P0", "red", "black"))
plot(cpm_counts[rownames(test_tissue)[1], ], col = ifelse(meta$tissue_type == "thymus", "red", "black"))
