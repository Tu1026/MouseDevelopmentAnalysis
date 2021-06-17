#This script will download curated metadata from google sheets.
#It then opens up fb_meta, an rds object saved by 01
#it then combines the two meta data frames

#Returns:
#The combined is saved as an rds object named complete_meta_data


library(tidyverse)
library(assertthat)
library(googlesheets4)

source(file = "Main/functions.R")

#-----Must change
#change sheet_id to your google sheet URL with the metadata on it.
sheet_id <- "1gFkSHD15wdd3FdCrZ51DQOi9RYQg01i1TXpIvkDVIz0"

#---download meta data containing replicate and dev stage info

replicate_meta <- read_sheet(sheet_id)

#---Open up fb_meta rds object

if (file.exists(paste0('Data/fb_meta.tsv'))) {
  fb_meta <- read_tsv(file = paste0('Data/fb_meta.tsv'))
} else {
  message("There was a problem reading fb_meta. If 01 was run successfully, it should be in Data/fb_meta")
}


#---join two tables

meta_data <- left_join(replicate_meta, fb_meta, by = c("id" = "File.accession"))
message("meta_data successfully downloaded")

#---save meta_data as tsv

write_tsv(x = meta_data, file = "Data/complete_meta_data.tsv")
