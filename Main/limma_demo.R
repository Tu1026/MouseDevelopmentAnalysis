## Demonstration of DE workflow using limma
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
## https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

# BiocManager::install("tximport")
# BiocManager::install("edgeR")

library(tidyverse)
library(limma)
library(edgeR)
library(tximport)

meta <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)

# count table files
files <- list.files("Data/Count_tables", ".tsv", full.names = TRUE)
names(files) <- meta$id

# use tximport pipeline to load all counts into a matrix

txi <- tximport(files, 
                type = "rsem",  # ENCODE pipeline uses RSEM
                txIn = FALSE,   # counts are already gene summarized
                txOut = FALSE,
                importer = read.delim)  # default read_tsv kept crashing

# only keep protein coding (must remove gene ID versions) 
count_mat <- txi$counts
rownames(count_mat) <- str_replace(rownames(txi$counts), "\\.[:digit:]+$", "")
count_mat <- count_mat[rownames(count_mat) %in% pc_mm$Gene_ID, ]

# design matrix for arbitrary comparison
group <- factor(meta$dev_stage == "P0")
design <- model.matrix(~group)

# remove genes that are not expressed using edgeR and scale normalization
# using TMM method (default)
dge <- DGEList(count_mat)
keep <- filterByExpr(dge, design)
dge <- dge[keep, ]
dge <- calcNormFactors(dge)

# convert counts to logCPM 
cpm_counts <- cpm(dge, log = TRUE, prior.count = 3)

# fit model with limma
fit <- lmFit(cpm_counts, design)
fit <- eBayes(fit, trend = TRUE)
de_table <- topTable(fit, number = Inf)

# include gene symbols and sort
pc_mm <- distinct(pc_mm, Symbol, .keep_all = TRUE)

de_table_clean <- de_table %>% 
  rownames_to_column(var = "Gene_ID") %>% 
  left_join(y = pc_mm[, c("Gene_ID", "Symbol")], by = "Gene_ID") %>% 
  arrange(desc(abs(logFC)))

# how many sig DE genes at FDR05?
sum(de_table_clean$adj.P.Val < 0.05)

# look at most different by logFC
max_fc <- de_table_clean$Gene_ID[which.max(de_table_clean$logFC)]
min_fc <- de_table_clean$Gene_ID[which.min(de_table_clean$logFC)]

plot(cpm_counts[max_fc, ])
plot(cpm_counts[min_fc, ])
