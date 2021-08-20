## This script ensures necessary libraries are installed, downloads list of 
## mouse protein coding genes, and saves the metadata file containing the 
## gene count file URLs.


# Package installation
# -----------------------------------------------------------------------------

packages <- c("tidyverse" "assertthat","googlesheets4", "skimr", "Hmisc", "corrplot", "pheatmap",
              "reshape2"
              )
install.packages(setdiff(packages, rownames(installed.packages())))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!"biomaRt" %in% rownames(installed.packages())) {
  BiocManager::install("biomaRt")
}

library(tidyverse)
library(assertthat)
library(biomaRt)
library(tools)
library(googlesheets4)


# Download metadata tsv containing gene count table URLs
# Taken from https://www.encodeproject.org/publication-data/ENCSR574CRQ/
# -----------------------------------------------------------------------------


url <- "https://www.encodeproject.org/documents/ab75e52f-64d9-4c39-aea0-15372479049d/@@download/attachment/ENCSR574CRQ_metadata.tsv"
file_path <- "Data/ENCSR574CRQ_metadata.tsv"

if (!file.exists(file_path)) {
  download.file(url, destfile = file_path)
}

stopifnot(is.readable(file_path))


# Get ensembl protein coding gene annotations
# -----------------------------------------------------------------------------


mm_mart <- "mmusculus_gene_ensembl"
mm_chromosome_filter <- c(1:19, "MT", "X", "Y")

# be aware of ensembl version! save and keep fixed if possible.

version <- str_extract(listMarts()$version[1], "[:digit:]+") 
mm_ensembl_outfile <- paste0("Data/ensembl_mouse_protein_coding_", version, ".tsv")

# get annotation tables, order to chromosomes, only take protein coding

mm_ensembl_mart <- useMart(biomart = "ensembl", dataset = mm_mart)

mm_attributes <- c(
  "chromosome_name",
  "transcription_start_site",
  "transcript_start",
  "transcript_end",
  "strand",
  "ensembl_gene_id",
  "mgi_symbol",
  "ensembl_transcript_id",
  "gene_biotype"
)

final_colnames <- c(
  "Chromosome",
  "Transcription_start_site",
  "Start",
  "End",
  "Strand",
  "Gene_ID",
  "Symbol",
  "Transcript_ID"
)

mm_anno_table <- getBM(
  attributes = mm_attributes,
  filter = "chromosome_name",
  values = mm_chromosome_filter,
  mart = mm_ensembl_mart,
  useCache = FALSE
)

# Only keeping protein coding genes

mm_protein_anno_table <- filter(mm_anno_table, gene_biotype == "protein_coding")
mm_protein_anno_table$gene_biotype <- NULL

# order the table by chromosome, then by TSS

mm_protein_anno_table <- mm_protein_anno_table %>%
  arrange(
    match(mm_protein_anno_table$chromosome_name, mm_chromosome_filter),
    transcription_start_site
  )

colnames(mm_protein_anno_table) <- final_colnames


# Save table
# -----------------------------------------------------------------------------


write.table(
  mm_protein_anno_table,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  mm_ensembl_outfile
)

stopifnot(is.readable(mm_ensembl_outfile))
