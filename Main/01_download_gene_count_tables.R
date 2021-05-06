## This script downloads gene count tables for the different mouse development
## transcriptome samples

library(tidyverse)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)
count_dir <- "Data/Count_tables"

if (!(dir.exists(count_dir))) {
  dir.create(count_dir)
}

# Focus on forebrain. Download gene count tables
# -----------------------------------------------------------------------------


fb_meta <- filter(file_meta,
             Biosample.term.name == "forebrain" &
             File.output.type == "gene quantifications")

# example on one file. implement iterative strategy.

id <- fb_meta$File.accession[1]
url <- fb_meta$S3.URL[1]
out_file <- paste0(count_dir, "/", id, ".tsv")

if (!file.exists(out_file)) {
  download.file(url, destfile = out_file)
}

# inspect

expr_dat <- read.delim(out_file, stringsAsFactors = FALSE)
gene_ids <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
table(gene_ids %in% pc$Gene_ID)
      