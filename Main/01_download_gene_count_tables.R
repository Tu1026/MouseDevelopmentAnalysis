## This script downloads gene count tables for the different mouse development
## transcriptome samples


library(tidyverse)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
file_meta <- read.delim("Data/ENCSR574CRQ_metadata.tsv", stringsAsFactors = FALSE)
count_dir <- "Data/Count_tables"

if (!(dir.exists(count_dir))) {
  warning("Creating Count_Dir")
  dir.create(count_dir)
}


#
# Focus on forebrain. Download gene count tables
# -----------------------------------------------------------------------------


fb_meta <- filter(file_meta,
             Biosample.term.name == "forebrain" &
             File.output.type == "gene quantifications")

# example on one file. implement iterative strategy.

dnld_data <- function(fb_meta)
  {
  for (i in 1:nrow(fb_meta)) {
    id <- fb_meta$File.accession[i]
    url <- fb_meta$S3.URL[i]
    out_file <- paste0(count_dir, "/", id, ".tsv")

    if (!file.exists(out_file)) {
      download.file(url, destfile = out_file)
    }
  }
}
fb_gen_ex <- list.files(count_dir)

if (is.null(fb_gen_ex)) {
  warning ("Expression data has not been downloaded yet. Downloading data...")
  dnld_data
}
#------md5sum

md5 <-lapply(X = paste0(count_dir,"/",list.files(count_dir)), FUN = md5sum)
names(md5)=c(fb_gen_ex)
md5_vector <- as.vector(unlist(md5,use.names=TRUE))
names(md5_vector) <- names(md5)
md5_vector

stopifnot(md5_vector == md5)
stopifnot(names(md5_vector) == names(md5))

names(md5_vector) <- str_replace(names(md5_vector),".tsv","")
md5_df <- data.frame(md5_vector)
md5_df["File.accession"] <- row.names(md5_df)
row.names(md5_df)=NULL
md5_df

dnld_md5 <- fb_meta%>%
  select(c(File.accession,md5sum))

md5_summary_table <- left_join(dnld_md5, md5_df, by = "File.accession")
colnames(md5_summary_table) = c("File.accession","downloaded_md5sum","generated_md5sum")

md5_summary_table <- md5_summary_table%>%
  mutate(equivalent = isTRUE(md5_summary_table$downloaded_md5sum==md5_summary_table$generated_md5sum))
# md5_summary_table

######HELP########## #TODO Troubleshoot why I'm getting False

stopifnot(md5_summary_table$Downloaded_md5sum == md5_summary_table$Generated_md5sum)


tryCatch(
  expr = {
    stopifnot(md5_summary_table$Downloaded_md5sum == md5_summary_table$Generated_md5sum)
  },
  error = function(cond){
    message("the md5 checksums do not match")
  }
)

  


#------iteratively get only pc genes with expression

# expr_dat <- read.delim(file = paste0(count_dir,"/","ENCFF227HKF.tsv"))
# gene_ids <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
# num_pc <- table(gene_ids %in% pc$Gene_ID)
# pc_gene_ex <- expr_dat %>%
#   select(c(gene_id,length,expected_count))%>%
#   mutate(gene_id = gene_ids) %>%
#   rename(Gene_ID = gene_id) %>%
#   filter(Gene_ID %in% pc$Gene_ID)%>%
#   left_join(y=pc, by="Gene_ID")%>%
#   select(c(Gene_ID,length,expected_count,Chromosome,Symbol))%>%
#   distinct()
# glimpse(pc_gene_ex)  


get_pc_genes <- function(gen_exp_files_names) 
  {
  if(!(dir.exists("Data/Processed_Count_tables"))){
    warning("Creating Data/Processed_Count_tables directory")
    dir.create("Data/Processed_Count_tables")
  }
    for (i in gen_exp_files_names) { 
      expr_dat <- read.delim(file = paste0(count_dir,"/",i))
      gene_ids <- str_replace(expr_dat$gene_id, "\\.[:digit:]", "")
      num_pc <- table(gene_ids %in% pc$Gene_ID)
      pc_gene_ex <- expr_dat %>%
        select(c(gene_id,length,expected_count))%>%
        mutate(gene_id = gene_ids) %>%
        rename(Gene_ID = gene_id) %>%
        filter(Gene_ID %in% pc$Gene_ID)%>%
        left_join(y=pc, by="Gene_ID")%>%
        select(c(Gene_ID,length,expected_count,Chromosome,Symbol))%>%
        distinct() %>%
        mutate(name = i)
      # stopifnot(gene_ids == pc_gene_ex$Gene_ID)
      # stopifnot(nrow(pc_gene_ex)==num_pc[2])
      # glimpse(pc_gene_ex)
      write_tsv(x=pc_gene_ex,
                file=paste0("Data","/","Processed_Count_tables","/","pc_gene_ex_",as.character(i)))
    }
}
my_expression_tables <- list.files(count_dir)

if (is.null("Data/Processed_Count_tables")){
  get_pc_genes(my_expression_tables)
}

#----------  open all processed count data

# open_pro_tables <- function() {
#   tables <- c()
#   for (i in my_expression_tables) {
#     print (i)
#     name <- paste0("pro_",i)
#     table <- read.delim(paste0("Data/Processed_Count_tables/pc_gene_ex_",i))
#     print(name)
#     glimpse(table)
#     name <-  table
#     append(x=tables, values=name)
#   return(tables)
#   }
# }
# 
# open_pro_tables()

pro_dir <- "Data/Processed_Count_tables"

processed_count_tables <- paste0(pro_dir,"/",list.files(pro_dir))
list_pro_count_tables <- lapply(processed_count_tables,read.delim)
names(list_pro_count_tables) <- list.files(pro_dir)

#todo make this more robust with regex
stopifnot (length(list_pro_count_tables)== length(my_expression_tables))

length(list_pro_count_tables)
nrows <-lapply(list_pro_count_tables, function(list_of_tables){nrow(list_of_tables)})
nrows <- unlist(nrows)
sum(nrows)
1*15015

nrow(list_pro_count_tables[[4]])

#---------- apply replicate metadata from google sheets

replicate_meta <- read_sheet("1gFkSHD15wdd3FdCrZ51DQOi9RYQg01i1TXpIvkDVIz0")



identify_in <- function(dataframe,replicate_meta){
  my_row <- replicate_meta[which(replicate_meta$file_dataset == dataframe$name[1]),]
  dataframe <- dataframe%>%
    mutate(tissue_type = factor(my_row[2]),
           dev_stage = factor(my_row[3]))
}

list_pro_count_tables <- lapply(list_pro_count_tables,identify_in,replicate_meta)
# list_pro_count_tables
# glimpse(list_pro_count_tables[1])
# head(list_pro_count_tables[1])


# replicate_meta[which(replicate_meta$file_dataset == "ENCFF227HKF.tsv"),]
# list_pro_count_tables[names(list_pro_count_tables)[1]]



#----------- conglomerate important info into 1 master table
create_master_table <- function(list_of_dataframes) {
  master_table <- data.frame()
  for (i in list_of_dataframes) {
    print(paste("Master_Table:", master_table))
    print(paste("list index:"),i)
    master_table <- rbind(master_table,i)
  return(master_table)
  }
}
master_table <- create_master_table(list_pro_count_tables)
nrow(master_table)

for (i in list_pro_count_tables){
  print (i)
}

master_table <- data.frame()
master_table <- rbind(master_table,list_pro_count_tables[[1]])
master_table <- rbind(master_table,list_pro_count_tables[[2]])

ncol(master_table)
nrow(master_table)
ncol(list_pro_count_tables[2])
stopifnot(nrow (master_table)==nrows)  

