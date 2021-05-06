library(tidyverse)

pc <- read.delim("~/Data/Metadata/ensembl_mouse_protein_coding_V98.tsv", stringsAsFactors = FALSE)
meta <- read.delim("~/scratch/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)
files <- read.delim("~/scratch/ENCSR574CRQ_files.txt")

fb <- filter(meta,
             Biosample.term.name == "forebrain" &
             File.output.type == "gene quantifications")

out_file <- "~/scratch/ENCFF895JXR.tsv"

if (!file.exists(out_file)) {
  download.file(fb[fb$File.accession == "ENCFF895JXR", "S3.URL"],
                destfile = out_file)
}

expr_dat <- read.delim(out_file, stringsAsFactors = FALSE)


tt <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
table(tt %in% pc$Gene_ID)

      