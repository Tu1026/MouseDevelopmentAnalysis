## This script downloads gene count tables for the different mouse development
## transcriptome samples. It also creates md5 checksum to confirm propper download

#---Returns: 
#Downloaded count tables,
# fb_meta,a forebrain filtered version of the total meta data


source(file = "Main/functions.R")

library(tools)
library(tidyverse)
library(assertthat)

count_dir <- "Data/Count_tables"
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)


#
#filter experimental meta_data for gene quantification, forebrain experiments
#

fb_meta <- filter(file_meta,
             Biosample.term.name == "forebrain" &
             File.output.type == "gene quantifications")

#---Save fb_meta as a rds
saveRDS(fb_meta, file = "Data/fb_meta")

#
#Creata a Data/Count_tables directory if one does not already exist. 
#
if (!(dir.exists(count_dir))) {
  warning("Creating Count_Dir")
  dir.create(count_dir)
}


#
#Deciding if dnld_data needs to be executed, or if data has already been downloaded
#

fb_gen_ex <- list.files(count_dir)

if (length(fb_gen_ex) != nrow(fb_meta)) {
    warning ("Expression data has not been downloaded yet. Downloading data...")
    dnld_data(fb_meta)
  } else { 
    warning ("Expression data has already been downloaded")
}



#
#Create a list of md5 checksums for each data table in l_expression_tables, and name according to count table names.
#


l_expression_tables_md5 <- lapply(paste0(count_dir,"/",list.files(count_dir)), create_md5)

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
warning("Md5 checksums were identical! Very nice")




 

