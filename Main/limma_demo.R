## Demonstration of DE workflow using limma
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
## https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
## https://bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html

# BiocManager::install("tximport")
# BiocManager::install("edgeR")
# BiocManager::install("limma")

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
                type = "rsem",  # ENCODE pipeline uses RSEM for gene counts
                txIn = FALSE,   # counts are already gene summarized
                txOut = FALSE,
                importer = read.delim)  # default read_tsv kept crashing

# only keep protein coding (must remove gene ID versions) 
count_mat <- txi$counts
rownames(count_mat) <- str_replace(rownames(txi$counts), "\\.[:digit:]+$", "")
count_mat <- count_mat[rownames(count_mat) %in% pc_mm$Gene_ID, ]
colnames(count_mat) <- paste(meta$id, meta$dev_stage, sep = "_")

# design matrix for dev stages
group <- factor(meta$dev_stage)
design <- model.matrix(~ 0 + group)


# remove genes that are not expressed using edgeR functionality
dge <- DGEList(count_mat, group = meta$dev_stage)
keep <- filterByExpr(dge, design)
dge <- dge[keep, ]

# gene counts of removed
removed <- rowSums(count_mat[!keep,])
count_mat[which.max(removed), ]

# scale normalization using TMM method (default)
dge <- calcNormFactors(dge)

# convert counts to log2CPM 
cpm_counts <- cpm(dge, log = TRUE, prior.count = 3)

# fit model with limma eBayes trend 
fit <- lmFit(cpm_counts, design)

# Define contrast - here, E10.5 versus all. Have to average contributions of 
# the comparison group

contr <- makeContrasts(
  E10_5 = groupE10.5 - (groupE11.5 + groupE12.5 + groupE13.5 + groupE14.5 + groupE15.5 + groupE16.5 + groupP0)/7,
  levels = design
)


contr <- makeContrasts(
  E10_5 = groupE10.5 - (groupE11.5 + groupE12.5 + groupE13.5 + groupE14.5 + groupE15.5 + groupE16.5 + groupP0)/7,
  levels = design
)


fit <- contrasts.fit(fit, contrasts=contr)
fit <- eBayes(fit, trend = TRUE)


# or, explicitly making a one versus all design matrix

group2 <- factor(meta$dev_stage == 'E10.5')
design2 <- model.matrix(~ group2)
fit2 <- eBayes(lmFit(cpm_counts, design2), trend = TRUE)

# look at top results
de_table <- topTable(fit, number = Inf)
de_table2 <- topTable(fit2, number = Inf)

# include gene symbols and sort
pc_mm <- distinct(pc_mm, Symbol, .keep_all = TRUE)

de_table_clean <- de_table2 %>% 
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
