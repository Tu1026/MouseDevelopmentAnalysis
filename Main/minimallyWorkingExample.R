#This entire script downloads meta data, encode data, and generate 3 flat 2-D matrices that have TPM of all samples

# Package installation
# -----------------------------------------------------------------------------

packages <- c("tidyverse", "assertthat","googlesheets4", "skimr", "Hmisc", "corrplot", "pheatmap",
              "reshape2", "sqldf"
)
install.packages(setdiff(packages, rownames(installed.packages())))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!"biomaRt" %in% rownames(installed.packages())) {
  BiocManager::install("biomaRt")
}

if (!"tximport" %in% rownames(installed.packages())) {
  BiocManager::install("tximport")
}

if (!"ENCODExplorer" %in% rownames(installed.packages())) {
  BiocManager::install("ENCODExplorer")
}

library(tidyverse)
library(assertthat)
library(biomaRt)
library(tools)
library(googlesheets4)
library(tximport)
source(file = "functions.R")


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



## This script downloads gene count tables for the different mouse development
## transcriptome samples. It also creates md5 checksum to confirm propper download

#---Returns: 
#Downloaded count tables,
# pre_meta,a forebrain filtered version of the total meta data





count_dir <- "Data/Count_tables"
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)



#filter experimental meta_data for gene quantification, forebrain experiments

# ini_set <- meta_data %>% filter(dev_stage == "E10.5") %>% select(tissue_type)
# 
# for (stage in uniq_dev_stages){
#   final_set <- ini_set %>% intersect(meta_data %>% filter(dev_stage == stage) %>% select(tissue_type)) 
# }
pre_meta <- filter(file_meta,
                    File.output.type == "gene quantifications")

#---Save pre_meta as a tsv

write_tsv(x = pre_meta, file = "Data/pre_meta.tsv")


#Creata a Data/Count_tables directory if one does not already exist. 

if (!(dir.exists(count_dir))) {
  warning("Creating Count_Dir")
  dir.create(count_dir)
}




gen_ex_fs <- list.files(count_dir)

if (length(gen_ex_fs) != nrow(pre_meta)) {
  message ("Expression data has not been downloaded yet. Downloading data...")
  dnld_data(pre_meta)
} else { 
  message ("Expression data has already been downloaded")
}




library(ENCODExplorer)


replicate_meta <- data.frame(id=character(), tissue_type=character(), dev_stage=character(), replicate=integer(), stringsAsFactors = FALSE)

for (idx in seq_along(unlist(pre_meta[,1]))){
   print(idx)
   experiment_meta = searchEncode(paste0(unlist(pre_meta[idx,1]), ' AND ', unlist(pre_meta[idx,2])))
   if (startsWith(unlist(experiment_meta[1,1]), "/experiments/")){
     idx_2 = 1
   }else {
     idx_2 = 2
   }
   tissue = experiment_meta$biosample_ontology.term_name[idx_2]
   print(tissue)
   dev_stage = experiment_meta$replicates.library.biosample.age[idx_2]
   print(dev_stage)
   if (dev_stage == 0){
     dev_stage = paste0("P", experiment_meta$replicates.library.biosample.age[idx_2])
   }else {
     dev_stage = paste0("E", experiment_meta$replicates.library.biosample.age[idx_2]) 
   }
   experiment_meta = searchEncode(unlist(pre_meta[idx,1]))
   replicate_num = strtoi(experiment_meta$biological_replicates[(which(experiment_meta$accession %in% unlist(pre_meta[idx,1])))])
   print(replicate_num)
   replicate_meta <- replicate_meta %>% add_row(id=unlist(pre_meta[idx,1]), tissue_type=tissue, dev_stage=dev_stage, replicate=replicate_num)
   
}

meta_data <- left_join(replicate_meta, pre_meta, by = c("id" = "File.accession"))

write_tsv(x = meta_data, file = "Data/complete_meta_data.tsv")




#---Open PC Table
pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv",
                 stringsAsFactors = FALSE)

#---Remove nondistinct Symbols and Gene_IDs from pc table
# and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  dplyr::select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)

write_tsv(x = pc_sub, file = "Data/pc_sub.tsv")



file_names_in_order <- paste(meta_data$id, ".tsv",sep="")
#---Read in all the count matrix
all_count_matrix <- tximport(paste0(count_dir,"/",file_names_in_order), type = "rsem", txIdCol = "gene_id", txIn = FALSE,
                        countsCol = "TPM", importer = read_delim)
colnames(all_count_matrix$abundance) <- file_names_in_order
#check the order of columns are correct
for (file in file_names_in_order){
  proper_path <- paste0(count_dir, "/", file)
  stopifnot(all(all_count_matrix$abundance[,file] == suppressMessages(read_tsv(proper_path))$TPM))
}
#-- Remove .tsv from tsv file name
colnames(all_count_matrix$abundance) <- str_replace(colnames(all_count_matrix$abundance), "\\.[:alpha:]*", "")
message("The columns are in the correct order in the imported matrix")

#colnames(all_count_matrix$abundance) <-  meta_data$id



rownames(all_count_matrix$abundance) <- str_replace(rownames(all_count_matrix$abundance),  "\\.[:digit:]*", "")
# Switch row names of the the count matrix to the gene symbol
old_row_names <- data.frame("Gene_ID" = rownames(all_count_matrix$abundance))
new_row_names <- old_row_names %>% left_join(pc_sub, by="Gene_ID") %>% filter(!is.na(Symbol))
all_count_matrix$abundance <- all_count_matrix$abundance[(rownames(all_count_matrix$abundance) %in% new_row_names$Gene_ID),]
rownames(all_count_matrix$abundance) <- new_row_names$Symbol


saveRDS(object = all_count_matrix$abundance, file = "Data/pc_count_matrix.rds")
write_tsv(x = data.frame(all_count_matrix$abundance), file = "Data/all_count_matrix.tsv")

odd_names <- colnames(all_count_matrix$abundance[,seq(1,ncol(all_count_matrix$abundance),2)])
even_names <- colnames(all_count_matrix$abundance[,seq(2,ncol(all_count_matrix$abundance),2)])
#---- Make sure that the pairs are indeed the same thing before doing pairwise operation
for (pair in 1:length(odd_names)){
  stopifnot((meta_data %>% filter(id == odd_names[pair]))["tissue_type"] == 
              (meta_data %>% filter(id == even_names[pair]))["tissue_type"])
}
for (pair in 1:length(odd_names)){
  stopifnot((meta_data %>% filter(id == odd_names[pair]))["dev_stage"] == 
              (meta_data %>% filter(id == even_names[pair]))["dev_stage"])
}
message("Pairs are the same good to go ahead")



avg_count_matrix <- (all_count_matrix$abundance[,seq(1,ncol(all_count_matrix$abundance),2)] + 
                       all_count_matrix$abundance[,seq(2,ncol(all_count_matrix$abundance),2)])/2;
saveRDS(object = avg_count_matrix, file = "Data/avg_pc_count_matrix.rds")

write_tsv(x = data.frame(avg_count_matrix), file = "Data/avg_count_matrix.tsv")

diff_count_matrix <- abs(all_count_matrix$abundance[,seq(1,ncol(all_count_matrix$abundance),2)] - 
                           all_count_matrix$abundance[,seq(2,ncol(all_count_matrix$abundance),2)])
saveRDS(object = diff_count_matrix, file = "Data/diff_pc_count_matrix.rds")
write_tsv(x = data.frame(diff_count_matrix), file = "Data/diff_count_matrix.tsv")
write_tsv(x = data.frame(rownames(diff_count_matrix)), "Data/rownames.tsv")

