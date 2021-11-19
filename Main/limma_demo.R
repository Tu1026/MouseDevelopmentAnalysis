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

# count table files - want to put in same order as meta for dev stages
files <- list.files("Data/Count_tables", ".tsv", full.names = TRUE)
file_ids <- str_replace(str_extract(files, "ENCFF.*"), ".tsv", "")
files <- files[order(match(file_ids, meta$id))]
names(files) <- meta$id

# use tximport pipeline to load all counts into a matrix

txi <- tximport(files, 
                type = "rsem",  # ENCODE pipeline uses RSEM for gene counts
                txIn = FALSE,   # counts are already gene summarized
                txOut = FALSE,
                importer = read.delim, # default read_tsv kept crashing
                countsCol = "expected_count")

# only keep protein coding (must remove gene ID versions) 
count_mat <- txi$counts
rownames(count_mat) <- str_replace(rownames(count_mat), "\\.[:digit:]+$", "")
count_mat <- count_mat[rownames(count_mat) %in% pc_mm$Gene_ID, ]
colnames(count_mat) <- paste(meta$id, meta$dev_stage, sep = "_")

# Convert gene ID rownames of count mat to gene symbols
pc_mm <- pc_mm %>% 
  filter(Symbol != "") %>% 
  distinct(Symbol, .keep_all = TRUE)

symbols <- left_join(data.frame(Gene_ID = rownames(count_mat)),
                     pc_mm[, c("Gene_ID", "Symbol")], by = "Gene_ID")

count_mat <- count_mat[symbols$Gene_ID, ]
stopifnot(identical(rownames(count_mat), symbols$Gene_ID))
rownames(count_mat) <- symbols$Symbol

# design matrix for dev stages
group <- factor(meta$dev_stage)
design <- model.matrix(~ 0 + group)
colnames(design) <- unique(meta$dev_stage)

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

# Define contrast of each dev stage vs all comparison. 
# Have to average contributions of the comparison group

contr <- makeContrasts(
  E10.5 = E10.5 - (E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E11.5 = E11.5 - (E10.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E12.5 = E12.5 - (E10.5 + E11.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E13.5 = E13.5 - (E10.5 + E11.5 + E12.5 + E14.5 + E15.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E14.5 = E14.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E15.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E15.5 = E15.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E16.5 + P0)/(n_distinct(meta$dev_stage)-1),
  E16.5 = E16.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + P0)/(n_distinct(meta$dev_stage)-1),
  P0 = P0 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5)/(n_distinct(meta$dev_stage)-1),
  levels = design
)


fit_contr <- contrasts.fit(fit, contrasts=contr)
fit <- eBayes(fit, trend = TRUE)
fit_contr <- eBayes(fit_contr, trend = TRUE)

# look at top results for each contrast

results <- lapply(colnames(design), function(x) {
  topTable(fit_contr, coef = x, number = Inf)
})
names(results) <- colnames(design)

topTable(fit_contr, coef = "E10.5", number = Inf)

# how many DE genes for each dev stage at FDR05
lapply(results, function(x) sum(x$adj.P.Val < 0.05))

hist(results$E10.5$P.Value, breaks = 100)

# This matrix shows the DE status of each gene for each contrast
contr_tests <- decideTests(fit_contr)

# example of genes that are never DE
no_de <- which(apply(contr_tests, 1, function(x) all(x == 0)))
contr_tests[no_de, ]

plot(cpm_counts["Tbx4", ])

# as expected by this design, no genes are DE in every contrast
all_de <- which(apply(contr_tests, 1, function(x) all(x == 1)))
contr_tests[all_de, ]

# arbitrary example of gene that is indicative of E10.5-E13.5. consider
# whether or not the genes must be non-DE/down-reg in the later time points
ix1 <- which(apply(contr_tests, 1, function(x) all(x[1:4] == 1)))
ix2 <- which(apply(contr_tests, 1, function(x) all(x[1:4] == 1 & x[5:8] %in% c(0, -1))))
contr_tests[ix1, ]
contr_tests[ix2, ]
plot(cpm_counts["Cdc45", ])

setdiff(ix1, ix2)
contr_tests[setdiff(ix1, ix2)[1], ]
plot(cpm_counts["Zfp276",])

# this ranks the most variable genes over all contrasts
topTable(fit_contr, number = 30)

plot(cpm_counts["Igf2bp1", ])
