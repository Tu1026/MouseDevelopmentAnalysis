#
#

source(file = "Main/functions.R")
source(file = "Main/libraries.R")
source(file = "Main/global_vars.R")
source("Main/01_download_gene_count_tables.R")
library(googlesheets4)

#This script filters meta_data and will download curated metadata from google sheets and combine it with fb_meta generated in 01. 
#The final meta_data table is named replicate_meta

#
#Filter for File.accession and add .tsv to each File.accession id
#
fb_meta_02 <- fb_meta %>%
  dplyr::select(File.accession) %>%
  transmute(File.accession = paste0(fb_meta$File.accession,".tsv"))
glimpse(fb_meta_02)

#
#download meta_data containing replicate and dev stage info
#

replicate_meta <- read_sheet("1gFkSHD15wdd3FdCrZ51DQOi9RYQg01i1TXpIvkDVIz0")

#
#join two tables
#
meta_data <- left_join(replicate_meta,fb_meta_02, by = "File.accession")
warning("meta_data successfully downloaded")

#
#Removal events
#


