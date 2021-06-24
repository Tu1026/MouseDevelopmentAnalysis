## This script does blah

library(tidyverse)

count_meta <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE)
count_dir <- "Data/Count_tables/"


open_count <- function(path) {
  
  # This function does blah
  
  read.delim(file = path,
             stringsAsFactors = FALSE,
             col.names = c(
               "Gene_ID",  # gene_id -> Gene_ID to be consistent with pcoding
               "transcript_id.s.",
               "length",
               "effective_length",
               "expected_count",
               "TPM",
               "FPKM",
               "posterior_mean_count",
               "posterior_standard_deviation_of_count",
               "pme_TPM",
               "pme_FPKM",
               "TPM_ci_lower_bound",
               "TPM_ci_upper_bound",
               "FPKM_ci_lower_bound",
               "FPKM_ci_upper_bound"
             ))
  
}


# Only want the forebrain samples

fb_ids <- count_meta[count_meta$tissue_type == "forebrain", "id"]
fb_files <- paste0(fb_ids, ".tsv")

stopifnot(all(fb_files %in% list.files(count_dir)))

# Load all forebrain count tables into a list

count_list <- lapply(paste0(count_dir, fb_files), open_count)
names(count_list) <- fb_ids

# Add ensembl protein coding symbols to each gene ID, if one exists. First
# must remove the version suffix on Gene_ID to make it consistent with the IDs
# in the pcoding table

proc_count_list <- get_symbols(count_list)

# Create a gene x sample matrix of counts

# using indexing from pre built mat
mat1 <- matrix(0, nrow = n_distinct(count_list$ENCFF302TQO$Gene_ID), ncol = length(count_list))
rownames(mat1) <- count_list$ENCFF302TQO$Gene_ID
colnames(mat1) <- fb_ids


# using do call on the list
mat2 <- do.call(cbind, lapply(count_list, `[[`, "TPM"))
rownames(mat2) <- count_list$ENCFF302TQO$Gene_ID
