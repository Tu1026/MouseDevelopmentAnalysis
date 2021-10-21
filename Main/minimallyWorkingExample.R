## This script ensures necessary libraries are installed, downloads list of 
## mouse protein coding genes, and saves the metadata file containing the 
## gene count file URLs.


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



## This script downloads gene count tables for the different mouse development
## transcriptome samples. It also creates md5 checksum to confirm propper download

#---Returns: 
#Downloaded count tables,
# fb_meta,a forebrain filtered version of the total meta data


source(file = "functions.R")

library(tools)
library(tidyverse)
library(assertthat)

count_dir <- "Data/Count_tables"
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)



#filter experimental meta_data for gene quantification, forebrain experiments

ini_set <- meta_data %>% filter(dev_stage == "E10.5") %>% select(tissue_type)

for (stage in uniq_dev_stages){
  final_set <- ini_set %>% intersect(meta_data %>% filter(dev_stage == stage) %>% select(tissue_type)) 
}
fb_meta <- filter(file_meta,
                    File.output.type == "gene quantifications")

#---Save fb_meta as a tsv

write_tsv(x = fb_meta, file = "Data/fb_meta.tsv")


#Creata a Data/Count_tables directory if one does not already exist. 

if (!(dir.exists(count_dir))) {
  warning("Creating Count_Dir")
  dir.create(count_dir)
}


#Deciding if dnld_data needs to be executed, or if data has already been downloaded+

browser()

fb_gen_ex <- list.files(count_dir)

if (length(fb_gen_ex) != nrow(fb_meta)) {
  message ("Expression data has not been downloaded yet. Downloading data...")
  dnld_data(fb_meta)
} else { 
  message ("Expression data has already been downloaded")
}



#Create a list of md5 checksums for each data table in l_expression_tables, and name according to count table names.


l_expression_tables_md5 <- lapply(paste0(count_dir,"/",list.files(count_dir)), md5sum)

count_dir_names <- list.files(count_dir)
names(l_expression_tables_md5) <- count_dir_names

#
#Checks. Check to ensure names are all in correct order. 
#

stopifnot(identical(count_dir_names, names(l_expression_tables_md5)))

#
#Select md5 from metadata. Create a new md5_meta table. Turn l_expression_tables_md5 into table and append to the md5_meta by File.accession
#

md5_meta <- fb_meta%>%
  dplyr::select(File.accession,md5sum) %>%
  mutate(File.accession = replace(File.accession, values = paste0(fb_meta$File.accession,'.tsv')))

t_expression_tables_md5 <- data.frame(File.accession = c(names(l_expression_tables_md5)),
                                      md5sum_downloaded = (as.vector(unlist(l_expression_tables_md5))
                                      )
)
stopifnot(all(md5_meta$File.accession %in% t_expression_tables_md5$File.accession))

md5_meta <- left_join(md5_meta,t_expression_tables_md5, "File.accession")

md5_meta <- md5_meta%>%
  mutate(md5sum_downloaded = as.character(unlist(md5sum_downloaded)))

#
#Check if md5 is identical
#

stopifnot(identical(md5_meta$md5sum, md5_meta$md5sum_downloaded))
message("Md5 checksums were identical! Very nice")

### Download and get replicate meta from BiocManager

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ENCODExplorer")

replicate_meta <- data.frame(id=character(), tissue_type=character(), dev_stage=character(), replicate=integer(), stringsAsFactors = FALSE)

for (idx in seq_along(unlist(fb_meta[,1]))){
   print(idx)
   experiment_meta = searchEncode(paste0(unlist(fb_meta[idx,1]), ' AND ', unlist(fb_meta[idx,2])))
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
   experiment_meta = searchEncode(unlist(fb_meta[idx,1]))
   replicate_num = strtoi(experiment_meta$biological_replicates[(which(experiment_meta$accession %in% unlist(fb_meta[idx,1])))])
   print(replicate_num)
   replicate_meta <- replicate_meta %>% add_row(id=unlist(fb_meta[idx,1]), tissue_type=tissue, dev_stage=dev_stage, replicate=replicate_num)
   
}

meta_data <- left_join(replicate_meta, fb_meta, by = c("id" = "File.accession"))

write_tsv(x = meta_data, file = "Data/complete_meta_data.tsv")





#This script does the following
#1) Open all of the count tables in Data/count_tables as a list named l_expr_tables
#2) Processes a single count table
#3) Generate Count-sample matrix

source('functions.R')

library(tidyverse)
library(assertthat)
library(stringr)
library(skimr)
library(Hmisc)
library(corrplot)
library(dplyr)


#---------------------------------------------------------------------------
########################## 1) Opening Count Tables
#--------------------------------------------------------------------------

# Opens all of the count tables in Data/count_tables as a list named l_expr_tables

#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

meta_data <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE, sep = "\t")

#---------------------------------------------------------------------------
# Global variables
#--------------------------------------------------------------------------
target_directory <- "Data/Count_tables"

#---------------------------------------------------------------------------
# Setting the order that the count tables are opened in
#--------------------------------------------------------------------------
# We are opening the count tables into a list. We want this list to be
# ordered in the same order as the meta_data table. This is to 
# represent the samples in order of develempental times

ordered_meta_data_id_names <- paste0(meta_data$id,".tsv")

#---------------------------------------------------------------------------
# Open the count tables
#--------------------------------------------------------------------------
# Open count tables and save inside a list. List should be ordered
# In the order of the meta_data table.
# Remove the .tsv after the names once the files are read

l_expr_tables <- lapply(ordered_meta_data_id_names,
                        open_expr_table,
                        target_directory)

names(l_expr_tables) <- str_replace(ordered_meta_data_id_names, ".tsv","")



#---------------------------------------------------------------------------
########################## 2) Generate  sample - expression data matrix
#--------------------------------------------------------------------------

# This script processes ALL count tables within l_expr_tables. 
# remove gene id versions from test data frames
# merge test data frame with pc table.
# Data frames are used to generate a  sample - expression data matrix
# sample - expression data matrix is saved


#---------------------------------------------------------------------------
# Open files
#--------------------------------------------------------------------------

# l_expr_tables is already opened: It is a list containing all of the count tables
# pc_sub is already generated: It is a processed pc table

#---Open meta_data
meta_data <- read.delim("Data/complete_meta_data.tsv",
                        stringsAsFactors = FALSE,
                        sep = "\t")

#---Open PC Table
pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv",
                 stringsAsFactors = FALSE)

#---Remove nondistinct Symbols and Gene_IDs from pc table
# and select only gene id and symbol, the relevant columns.
pc_sub <- pc %>% 
  dplyr::select(Gene_ID, Symbol)%>%
  distinct(Gene_ID, .keep_all = TRUE) %>% 
  distinct(Symbol, .keep_all = TRUE)

#---------------------------------------------------------------------------
# Prep PC and l_expr_tables for merging with each other, then merge
#--------------------------------------------------------------------------

#---Run remove_gene_id_vers on count tables to remove gene id versions
stopifnot("Gene_ID" %in% colnames(l_expr_tables[[1]]))
l_expr_tables_noid <- lapply(l_expr_tables, remove_gene_id_vers)

#---merge pc_sub table with all count tables
l_expr_merged <- lapply(l_expr_tables_noid, merge_count_pc, pc_sub)

#---------------------------------------------------------------------------
# Get only Protein Coding genes from merged expression - pc data frames. 
#---------------------------------------------------------------------------
# Now that our count tables have Protein Symbols, we can 
# Use the presence of a symbol to filter for protein coding genes 
l_expr_symbols <- lapply(l_expr_merged, get_genes_with_symbols)

#---------------------------------------------------------------------------
# Creating the sample - expression data matrix
#--------------------------------------------------------------------------
#---Test to see that each frame in l_expr_symbols has the same order of symbols
# If this is true, then we can just iteratively add the lists onto a matrix

stopifnot(test_orders(l_expr_symbols))


#---Create empty expression matrix and fill with expression data
count_matrix <- matrix(data = 0,
                       nrow = nrow(l_expr_symbols$ENCFF227HKF),
                       ncol = length(l_expr_symbols))

rownames(count_matrix) <- l_expr_symbols$ENCFF227HKF$Symbol
colnames(count_matrix) <- names(l_expr_symbols)

for (colname in colnames(count_matrix)) { 
  count_matrix[ , colname] <- l_expr_symbols[[colname]][["TPM"]]
}

#---Save expression matrix as rds

saveRDS(object = count_matrix, file = "Data/pc_count_matrix.rds")

#---------------------------------------------------------------------------
# Building an AVG Matrix
#---------------------------------------------------------------------------
# As seen in meta_data, every mouse developmental stage has 2 samples assosiated
# with it. We want to create a new matrix with those replicates averaged so
# we can easily perform analysis on how our counts change over the dev stages
# Build a new matrix containing the averages of counts across replicates
# Called avg_count_matrix


## 30 12 16 16 16 24 24 18 number of samples for each dev stage


## Collapses the sampels into average
avg_count_matrix <- (count_matrix[,seq(1,ncol(count_matrix),2)] + count_matrix[,seq(2,ncol(count_matrix),2)])/2;

# avg_count_matrix <- matrix(data = 0,
#                            nrow = nrow(count_matrix),
#                            ncol = ncol(count_matrix)/6)
rownames(avg_count_matrix) <- rownames(count_matrix)

list_of_avg_matrices <- list()

for (tissue in unique(meta_data$tissue_type)){
  generic_mat <- matrix(nrow=nrow(avg_count_matrix), ncol = lengths(meta_data%>%filter(tissue_type==tissue))[1]/2)
  walk_count <- 1
  dev_stages <- c()
   for (exp_id in colnames(avg_count_matrix)){
     if ((meta_data %>% filter(id==exp_id))$tissue_type == tissue){
       dev_stages <- c(dev_stages, (meta_data %>% filter(id==exp_id))$dev_stage)
       generic_mat[,walk_count] <- avg_count_matrix[,exp_id]
       walk_count <- walk_count + 1
     }
   }
         colnames(generic_mat) <- dev_stages
         rownames(generic_mat) <- rownames(avg_count_matrix)

       list_of_avg_matrices[[tissue]] <- generic_mat
}





saveRDS(object = list_of_matrices, file = "Data/avg_pc_count_matrix.rds")


#---------------------------------------------------------------------------
# Creating a difference count matrix
#---------------------------------------------------------------------------
# We want to identify which genes in a pair of replicates differ the most
# in count value.

## 30 12 16 16 16 24 24 18 number of samples for each dev stage




list_of_diff_matrices <- list()
for (tissue in unique(meta_data$tissue_type)){
  generic_mat <- matrix(nrow=nrow(count_matrix), ncol = lengths(meta_data%>%filter(tissue_type==tissue))[1]/2)
  walk_count <- 1
  dev_stages <- c()
  for (mydev_stage in uniq_dev_stages){
    meta_subset <- filter(meta_data, meta_data$dev_stage == mydev_stage, tissue_type==tissue)
    if (nrow(meta_subset) != 0){
      dev_stages <- c(dev_stages, mydev_stage)
      generic_mat[,walk_count] <- abs(count_matrix[ , meta_subset$id[1]] - count_matrix[ , meta_subset$id[2]])
      walk_count <- walk_count + 1
    }
  }
  colnames(generic_mat) <- dev_stages
  rownames(generic_mat) <- rownames(avg_count_matrix)
  list_of_diff_matrices[[tissue]] <- generic_mat
}



#---Save avg expression matrix as rds

saveRDS(object = diff_count_matrix, file = "Data/diff_pc_count_matrix.rds")
