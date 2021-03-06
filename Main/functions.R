merge_tables <- function(pc_table,expression_df) {
  #
  #joins the protein coding table onto the expression coding table by gene_id = Gene_Id.
  #
  #@param: pc_table: the protein coding table downlaoded by ensemble. Should be named "pc"
  #@param: expression_df: A given expression data frame. Should be within "l_expr_tables"
  #return: new data frame thats merged both tables merged by "gene_id" columns
  
  merged_table <- left_join(expression_df , pc_table, by = c("Gene_Id" = "gene_id"))
  return (merged_table)
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

test_orders <- function(l_expr_symbols) { 
  reference_frame <- l_expr_symbols[[1]]$Symbol 
  for (frame_index in 1:length(l_expr_symbols)) {
    test_frame <- l_expr_symbols[[frame_index]]$Symbol
    stopifnot(all(reference_frame == test_frame))
  }
  message("The order of all symbols within l_expr_symbols is the same")
  return(TRUE)
}

# #
# #Create md5sum for each downloaded expression table
# #
# 
# create_md5 <- function(path_to_dataset) {
#   #
#   #Creates an md5check sum for an inputted dataset.
#   #
#   #@param: In string, the path to the datset you want an md5checksum for
#   #     ex: "Data/Count_tables/ENCFF9760LT.tsv"
#   #@return: an md5checksum for inputted dataset
#   
#   md5 <- md5sum(path_to_dataset)
#   return(md5)
# }


#
#lapply to create a list of all expression tables. l_expr_tables is a list containing all expression tables
#

create_names_of_tables<- function(count_dir) {
  #
  #Create a list of all the names of the tables within count_dir
  #
  #@param: directory where expression tables are
  #@return: list of all the names of your expression tables. Includes location in directory
  my_list <- list()
  for (i in list.files(count_dir)) { 
    my_list <- append(my_list, paste0(count_dir,"/",i))
  }
  return (my_list)
}


remove_gene_id_vers <- function(data_frame) {
  #
  #Removes the version ids of the gene_id column within all expr tables
  #
  #@param : data frame for expression
  #@return : data frame with gene_id versions removes
  data_frame$Gene_ID <- str_replace(data_frame$Gene_ID, "\\.[:digit:]*", "")
  return(data_frame)
}




#
#Generating Positive and Negative Checks for remove_gene_id_vers()
#
check_gene_id_vers_removal <- function(dataframe) {
  #
  #Checks to see that remove_gene_id_vers was successful
  #
  #@param: Dataframe
  #@return: TRUE= Removal was successfull. FALSE= removal failed
  my_checks <- c( "ENSMUSG00000021028.7  ENSMUSG00000032667.13 ENSMUSG00000016758.3  ENSMUSG00000014418.11 ENSMUSG00000092505.1  ENSMUSG00000045395.3 
                ENSMUSG00000084007.3  ENSMUSG00000087052.1  ENSMUSG00000096615.1  ENSMUSG00000081432.1  ENSMUSG00000076888.1  ENSMUSG00000071816.5 
                ENSMUSG00000064860.1  ENSMUSG00000062098.7  ENSMUSG00000093867.3  ENSMUSG00000027257.9  ENSMUSG00000093806.1  ENSMUSG00000087496.1 
                ENSMUSG00000089999.1  ENSMUSG00000029072.3  ENSMUSG00000101458.1  ENSMUSG00000025747.8  ENSMUSG00000095015.3  ENSMUSG00000090383.1 
                ENSMUSG00000057461.1  ENSMUSG00000081266.1  ENSMUSG00000076485.2")
  my_negative_checks <- str_extract_all(my_checks, boundary("word")) #After running remove_gene_id_vers on a data table, this check should yield false
  my_positive_checks <- str_replace(my_negative_checks[[1]], "\\.[:digit:]*", "") #After running remove_gene_id_vers on a data table, this check should yield TRUE
  if (all(!my_negative_checks %in% dataframe$gene_id) == TRUE & all(my_positive_checks %in% dataframe$gene_id) == TRUE) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}


#This script takes l_pc_expr and generates a bunch of summary stat data
generate_summary_stats <- function(data_frame) {
  #
  #Purpose: Take a signle data frame of l_pc_expr and generate a bunch of summary information and processed dataframes
  #
  #@param: data_frame: A df of l_pc_expr
  #@return: A list containing summary information, and data frames
  #
  #---Filter out transcripts that do have a protein symbol
  df_sym <- data_frame %>%
    filter(!(is.na(data_frame$Symbol)))
  #---Check that all of your transcriopts do indeed have symbols
  stopifnot(!(all(is.na(df_sym$Symbol))))
  
  #---Filter out transcripts that NOT have a protein symbol
  df_nosym <- data_frame %>%
    filter((is.na(data_frame$Symbol)))
  #---Check that all of your transcriopts do  NOT have symbols
  stopifnot((all(is.na(df_nosym$Symbol))))
  
  #Check to make sure your number of symbol transcripts + non-symbol transcripts adds up to total
  stopifnot(nrow(df_sym)+nrow(df_nosym) == nrow(data_frame))
  
  #---Filter For the protein coding transcripts with expected count 0
  df_sym_ec0 <- df_sym %>%
    filter(df_sym$expected_count == 0)
  #---Filter for Protein coding transcripts that have expected count >0
  df_sym_ec1 <- df_sym %>%
    filter(df_sym$expected_count > 0)
  #---Check to make sure the number of transcripts with expression 0 + number of transcripts with expression >0 is equal to the total number of PC transcripts
  stopifnot(nrow(df_sym_ec0)+nrow(df_sym_ec1)==nrow(df_sym))
  
  #---Filter for the number of unique pc transcripts in total (expressed + non expressed)
  n_genes <- df_sym%>%
    group_by(Symbol) %>%
    tally()%>%
    nrow()
  #---For Expressed PC Transcripts. Get the number of Symbols (Unique Genes)
  n_genes_expr <- df_sym_ec1%>%
    group_by(Symbol) %>%
    tally()%>%
    nrow()
  #---Filter for the number of unique genes NOT expressed
  n_genes_not_expr <- df_sym_ec0%>%
    group_by(Symbol) %>%
    tally()%>%
    nrow()
  #---Check to make sure the number of unique genes expressed + number of unique genes not expressed == number of unique transcripts with Symbols
  stopifnot(n_genes_expr + n_genes_not_expr == n_genes)
  
  current_list <- list(data_frame, df_sym, df_nosym, df_sym_ec0, df_sym_ec1, n_genes, n_genes_expr, n_genes_not_expr)
  names(current_list) <- c("data_frame", "df_sym", "df_nosym", "df_sym_ec0", "df_sym_ec1", "n_genes", "n_genes_expr", "n_genes_not_expr")
  return(current_list)
}


generate_filter_plot <- function(list_of_sum_stats) { 
  #
  #list_of_sum_stats is generated in 09. This script will generate plots for that information. 
  #It is about the merged pc_expression tables and the fitering that happens on it
  #
  
  #names should be a character vector of all the names inside our list
  names <- c("Merged Transcripts","Transcripts with Symbol(PC Genes)","Transcripts with no Symbol" , "Non-Expressed PC Transcripts", "Expressed PC Transcripts",
             "Number of Genes in PC table", "Expressed Genes", "Non Expressed Genes")
  values <- c(nrow(list_of_sum_stats$data_frame), nrow(list_of_sum_stats$df_sym), nrow(list_of_sum_stats$df_nosym), nrow(list_of_sum_stats$df_sym_ec0),
              nrow(list_of_sum_stats$df_sym_ec1), list_of_sum_stats$n_genes, list_of_sum_stats$n_genes_expr, list_of_sum_stats$n_genes_not_expr)
  df <- data.frame(names = factor(names, levels = names), 
                   values = values)
  return (ggplot(df) + geom_col(mapping = aes(x=names, y=values)))
}


open_expr_table <- function(expression_table.tsv,
                            target_directory,
                            opening_order) {
  #open a single  experession table. 
  #expression_table.tsv = the name of the expression table you wish to open as string
  #target_directory = the directory where the expression table is saved, as string
  #column_names = the column names you wish to open the expression table as. 
  #Setting this as a vector will ovverwite the default names, a vector strings
  #opening_order is a 
  
  table <- read.delim(file = paste0(target_directory,"/",expression_table.tsv), 
                          sep = "\t",
                          col.names = c(
                            "Gene_ID",  # gene_id -> Gene_ID to be consistent with pcoding
                            "transcript_id.s.",
                            "length",
                            "effective_length",
                            "expected_count",
                            "TPM",
                            "FPKM",
                            "posterior_mean_count",
                            "posterior_standard_deviation_of_count",
                            "pme_TPM",
                            "pme_FPKM",
                            "TPM_ci_lower_bound",
                            "TPM_ci_upper_bound",
                            "FPKM_ci_lower_bound",
                            "FPKM_ci_upper_bound"))
      return(table)
}

merge_count_pc <- function(count_table, pc_sub) { 
  #merge a count table onto the pc table
  merged_table <- count_table %>%
    left_join(y = pc_sub, by = "Gene_ID") %>% 
    relocate(Symbol, .after = Gene_ID)
  return (merged_table)
}

get_genes_with_symbols <- function(merged_table) {
  #---Return a table containing the transcripts that have a protein symbol
  with_symbols <- filter(merged_table,!(is.na(Symbol)))
  return (with_symbols)
}


rowSD <- function(matrix) {
  #
  # For a matrix, return the st_dev of each row
  #
  sd_matrix <- matrix(0,
                      nrow = nrow(matrix),
                      ncol = 1)
  rownames(sd_matrix) <- rownames(matrix)
  
  for (row_num in 1:nrow(matrix)) {
    row <- matrix[row_num,]
    st_dev <- sd(row)
    sd_matrix[row_num,] <- st_dev
  }
  return(sd_matrix[,1])
}

rowVar <- function(matrix) {
  #
  # For a matrix, return the variance of each row
  #
  var_matrix <- matrix(0,
                      nrow = nrow(matrix),
                      ncol = 1)
  rownames(var_matrix) <- rownames(matrix)
  
  for (row_num in 1:nrow(matrix)) {
    row <- matrix[row_num,]
    row_var<- var(row)
    var_matrix[row_num,] <- row_var
  }
  return(var_matrix[,1])
}

