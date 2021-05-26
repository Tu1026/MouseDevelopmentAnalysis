## This script downloads gene count tables for the different mouse development
## transcriptome samples

whats_run <- c(whats_run,"01")

library(tidyverse)
library(tools)
library(googlesheets4)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)
count_dir <- "Data/Count_tables"


#
# Focus on forebrain. Download gene count tables
# -----------------------------------------------------------------------------


#
#filter experimental meta_data for gene quantification, forebrain experiments
#
fb_meta <- filter(file_meta,
             Biosample.term.name == "forebrain" &
             File.output.type == "gene quantifications")
glimpse(fb_meta)

#
#Creata a Data/Count_tables directory if one does not already exist. 
#
if (!(dir.exists(count_dir))) {
  warning("Creating Count_Dir")
  dir.create(count_dir)
}


dnld_data <- function(fb_meta) {
  #
  #Downloads expression data by using the URL provided in fb_meta
  #Only downloads gene quantification data for forebrain experiments
  #
  #@param fb_meta: The experimental meta_data containing only quantification data for forebrain experiments
  #@return: Downloads expression data into the count_dir directory (Cata/Count_tables)
  
  for (i in 1:nrow(fb_meta)) {
    id <- fb_meta$File.accession[i]
    url <- fb_meta$S3.URL[i]
    out_file <- paste0(count_dir, "/", id, ".tsv")

    if (!file.exists(out_file)) {
      download.file(url, destfile = out_file)
    }
  }
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
#Create md5sum for each downloaded expression table
#

create_md5 <- function(path_to_dataset) {
  #
  #Creates an md5check sum for an inputted dataset.
  #
  #@param: In string, the path to the datset you want an md5checksum for
  #     ex: "Data/Count_tables/ENCFF9760LT.tsv"
  #@return: an md5checksum for inputted dataset
  
  md5 <- md5sum(path_to_dataset)
  return(md5)
}



#
#Create a list of md5 checksums for each data table in l_expression_tables, and name according to count table names.
#


l_expression_tables_md5 <- lapply(paste0(count_dir,"/",list.files(count_dir)),create_md5)

count_dir_names <- list.files(count_dir)
names(l_expression_tables_md5) <- count_dir_names

#
#Checks. Check to ensure names are all in correct order. 
#

stopifnot(identical(count_dir_names, names(l_expression_tables)))
stopifnot(identical(count_dir_names, names(l_expression_tables_md5)))

#
#Select md5 from metadata. Create a new md5_meta table. Turn l_expression_tables_md5 into table and append to the md5_meta by File.accession
#

#TODO: Why does as.list() work here but not list(). Lists in R are complicated thats for sure

md5_meta <- fb_meta%>%
  dplyr::select(File.accession,md5sum) %>%
  mutate(File.accession = replace(File.accession, values = paste0(fb_meta$File.accession,'.tsv')))
  
glimpse(md5_meta)  

t_expression_tables_md5 <- data.frame(File.accession = c(names(l_expression_tables_md5)),
                                      md5sum_downloaded = (as.vector(unlist(l_expression_tables_md5))
                                      )
)

glimpse(t_expression_tables_md5)

stopifnot(all(md5_meta$File.accession %in% t_expression_tables_md5$File.accession))

md5_meta <- left_join(md5_meta,t_expression_tables_md5, "File.accession")

glimpse(md5_meta)

md5_meta <- md5_meta%>%
  mutate(md5sum_downloaded = as.character(unlist(md5sum_downloaded)))

glimpse(md5_meta)

#
#Check if md5 is identical
#

stopifnot(identical(md5_meta$md5sum, md5_meta$md5sum_downloaded))
warning("Md5 checksums were identical! Very nice")




# pc1 <- pc[1:10,]
# pc2 <- pc[1:10,]
# pc1$Chromosome == pc2$Chromosome
# all(pc1$Chromosome == pc2$Chromosome)
# identical(pc1$Chromosome, pc2$Chromosome)

#
#Removal events
#

rm("create_md5","dnld_data","t_expression_tables_md5","md5_meta","l_expression_tables_md5","fb_gen_ex")



 

