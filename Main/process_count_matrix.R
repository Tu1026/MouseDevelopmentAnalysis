library(tidyverse)
library(assertthat)
library(stringr)
library(Hmisc)
library(corrplot)


#---------------------------------------------------------------------------
# Open Correlation Matrix and meta_data
#---------------------------------------------------------------------------

count_matrix <- readRDS(file = "Data/pc_count_matrix.rds")
avg_count_matrix <- readRDS(file = "Data/avg_pc_count_matrix.rds")
diff_count_matrix <- readRDS(file = "Data/diff_pc_count_matrix.rds")


#---------------------------------------------------------------------------
# Generating Pearson Correlation Matrix
#---------------------------------------------------------------------------
# Generating pearson correlation matrix to see correlation between sample's gene
# counts

# Also generate correlation matrix's p-values

pear_correlation_matrix_rcorr_object <- rcorr(count_matrix, type = "pearson")
pear_correlation_matrix <- pear_correlation_matrix_rcorr_object$r
pear_correlation_matrix_pvalues <- pear_correlation_matrix_rcorr_object$P

#---Generating heatmap for sample correlation matrix

pear_heatmap <- heatmap(x = pear_correlation_matrix,
                   sym = TRUE,
                   Rowv = NA,
                   Colv = NA, 
                   revC= TRUE)

#---------------------------------------------------------------------------
# Generating spearman Correlation Matrix
#---------------------------------------------------------------------------
# Generating spearman correlation matrix to see correlation between sample's gene
# counts if they were ordered

# Also generate correlation matrix's p-values

spear_correlation_matrix_rcorr_object <- rcorr(count_matrix, type = "spearman")
spear_correlation_matrix <- spear_correlation_matrix_rcorr_object$r
spear_correlation_matrix_pvalues <- spear_correlation_matrix_rcorr_object$P

#---Generating heatmap for sample correlation matrix

spear_heatmap <- heatmap(x = spear_correlation_matrix,
                        sym = TRUE,
                        Rowv = NA,
                        Colv = NA, 
                        revC= TRUE)


#---Check to confirm the averaging functioned properly

for (check_number in 1:5) {
  check_name <- sample(rownames(avg_count_matrix), size = 1)
  
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,1:2]), 
                                  avg_count_matrix[check_name, 1]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,3:4]), 
                                  avg_count_matrix[check_name, 2]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,9:10]),
                                  avg_count_matrix[check_name, 5]))
  stopifnot(assertthat::are_equal(mean(count_matrix[check_name,15:16]), 
                                  avg_count_matrix[check_name, 8]))
}

#---------------------------------------------------------------------------
# Comments on Spearman and Pearson Correlation Matrices
#---------------------------------------------------------------------------

# The pearson correlation decreases as the mouse samples progress through
# the dev stages

# On the other hand, the spearman correlation stays relatively constant
# as the mouse samples progress through the dev stages

# This is likely due to the fact that the spearman correlation is correlating
# pair-wise values which are paired with respect to count.
# Whereas the pearson is correlating pair-wise values which are paired with 
# respect to the Gene Symbols

# These observations suggest that the expression of several genes changes
# throughout the developmental stages. However, the magnitude of expression
# remains relatively constant


#---------------------------------------------------------------------------
# Comparing the count values of two replicate samples
#---------------------------------------------------------------------------

# Generate a scatterplot that compares the counts of two replicate samples

Replicate_1 <- count_matrix[,1]
Replicate_2 <- count_matrix[,2]

replicate_df <- data.frame(Replicate_1 = Replicate_1,
                           Replicate_2 = Replicate_2, 
                           row.names = rownames(Replicate_1))

replicate_scatterplot <- ggplot(replicate_df)+ 
  geom_point(mapping = aes(x = Replicate_1, y = Replicate_2))+
  labs( title = "Pair-wise counts of two replicate samples")

replicate_scatterplot


#---------------------------------------------------------------------------
# Comparing the count values of two non-replicate samples
#---------------------------------------------------------------------------

# Generate a scatterplot that compares the counts of two replicate samples

non_Replicate_1 <- count_matrix[,1]
non_Replicate_2 <- count_matrix[,16]

non_replicate_df <- data.frame(non_Replicate_1 = non_Replicate_1,
                               non_Replicate_2 = non_Replicate_2, 
                                row.names = rownames(Replicate_1))

non_replicate_scatterplot <- ggplot(non_replicate_df)+ 
  geom_point(mapping = aes(x = non_Replicate_1, y = non_Replicate_2))+
  labs( title = "Pair-wise counts of two non replicate samples")

non_replicate_scatterplot


#---------------------------------------------------------------------------
# Comments on replicate scatterplot and non-replicate scatterplot
#---------------------------------------------------------------------------

# As expected, the scatterplot of the replicates showed a strong correlation.
# Whereas the non-replicate scatterplot was much more scattered,
# indicating a weaker correlation
# This is indicative of the calculated pearson correlation coefficients,
# which were 0.99 and 0.60 for the two graphs respectively

#---------------------------------------------------------------------------
# Identifying differentially expressed genes between replicates
#---------------------------------------------------------------------------

# Ideally, replicates should have similar counts for all genes. 
# Here we will explore the genes that dont and add them to alist

# --- Looking at E105 first
E105 <- diff_count_matrix[, colnames(diff_count_matrix)[1]]

E105_df <- data.frame("E105" = E105)
E105_df <- arrange(E105_df, desc(E105_df$E105))

E105_hist <- ggplot(E105_df, aes(x = E105)) + 
  geom_histogram(binwidth = 10) +
  labs( title = "E105 most differentially expresed genes") +
  scale_y_continuous(trans = "log10")
E105_hist

l_most_diff_genes <- list("E105" = head(E105_df, n = 10))

# --- Doing the same for all dev stages
l_most_diff_genes <- list()


get_most_diff <- function(column) {
  my_df <- data.frame("Diff_count" = column)
  arranged_df <- arrange(my_df, desc(my_df))
  return(top_n(arranged_df, n = 10))
}

most_diff_genes <- apply(diff_count_matrix, MARGIN = 2, FUN = get_most_diff)


#---------------------------------------------------------------------------
# Identifying most expressed genes across samples
#---------------------------------------------------------------------------

# Here we will identify the genes that are most expressed in each sample

# --- Looking at E105 first
E105 <- avg_count_matrix[, colnames(avg_count_matrix)[1]]

E105_df <- data.frame("E105" = E105)
E105_df <- arrange(E105_df, desc(E105_df$E105))

l_most_expr_genes <- list("E105" = head(E105_df, n = 10))

# --- Doing the same for all dev stages
l_most_expr_genes <- list()


get_most_expr<- function(column) {
  my_df <- data.frame("count" = column)
  arranged_df <- arrange(my_df, desc(my_df))
  return(top_n(arranged_df, n = 10))
}

most_expr_genes <- apply(avg_count_matrix, MARGIN = 2, FUN = get_most_expr)


#---------------------------------------------------------------------------
# Comments on most expr and most diff genes
#---------------------------------------------------------------------------

# It looks like Tuba1a is the most expressed gene
# The problem is that many of the genes that are highly expressed
# are also seen to have the most deviation between replicates
# so we have to ask ourselves how reliable these counts are


#---------------------------------------------------------------------------
# Exploring genes that are never expressed
#---------------------------------------------------------------------------

# THere are some genes that never got counted across all dev stages
# Here I've arbitrarily chosen <=1 as never expressed. May chance later


never_expressed_genes <- avg_count_matrix[  avg_count_matrix[ , 1] <= 1 &
                        avg_count_matrix[ , 2] <= 1 & 
                        avg_count_matrix[ , 3] <= 1 &
                        avg_count_matrix[ , 4] <= 1 &
                        avg_count_matrix[ , 5] <= 1 &
                        avg_count_matrix[ , 6] <= 1 &
                        avg_count_matrix[ , 7] <= 1 &
                        avg_count_matrix[ , 8] <= 1 , ]

#---------------------------------------------------------------------------
# Defining a lowly expressed gene
#---------------------------------------------------------------------------

# Maybe a lowly expressed gene is in the 10th percentile for expression


# It could also be 2 stdeviations below the mean of expression

avg_count_matrix 






# lowly_expressed_genes <- matrix(,
#                                 nrow = 10000,
#                                 ncol = ncol(avg_count_matrix)
# )
# 
# 
# colnames(lowly_expressed_genes) <- colnames(avg_count_matrix)
                                
# for (column_name in colnames(lowly_expressed_genes)) {
#   lowly_expressed_genes[, column_name] <- 
#     avg_count_matrix[avg_count_matrix[, column_name] <= 1, column_name]
# }



#---------------------------------------------------------------------------
# Compressing Metadata 
#---------------------------------------------------------------------------
# Not sure if there are betters ways to extract specific columns from matrices,
# But this is the best I could do
# I'm literally just trying to figure out how to link the replicates onto the expression matrix columns

matrix_c_names <- colnames(count_matrix)


meta_data_group <- meta_data %>%
  select(id, dev_stage) %>%
  group_by(dev_stage)%>%
  dplyr::summarize(ids = str_split(paste(id, collapse = ","), pattern = ","))
rownames(meta_data_group) <- meta_data_group$dev_stage

indexes <- list()
for (i in meta_data_group$dev_stage) {
  row <- meta_data_group[i,]
  names <- unlist(row$ids)
  cur_indexes <- which(matrix_c_names == names)
  indexes <- append(indexes,list(cur_indexes))
}

indexes_almost_merged <- do.call(rbind, indexes)
meta_data_group <- meta_data_group %>%
  mutate(index1 = indexes_almost_merged[,1],
         index2 = indexes_almost_merged[,2])

#---------------------------------------------------------------------------
# Comparing reads between replicates
#---------------------------------------------------------------------------


#-------doing for 1

e105 <- meta_data_group["e10.5",]
e105 <- unlist(e105$ids)
indexes <- which(matrix_col_names == e105)

extract_matrix_columns <- function(index, count_matrix) {
  return (count_matrix[,index])
}

e105_count_matrix <- do.call(cbind, lapply(indexes, 
                                           extract_matrix_columns,
                                           count_matrix))

replicate_differences <- abs(e105_count_matrix[,1] - e105_count_matrix[,2])

#-------Doing for all count tables


master_diff_matrix <- matrix(data = 0,
                             nrow = nrow(count_matrix), 
                             ncol = nrow(meta_data_group),   
                             dimnames =  list(rownames(count_matrix), meta_data_group$dev_stage))
for (i in 1:nrow(meta_data_group)) {
  row <- meta_data_group[i,]
  matrixcol1 <- count_matrix[,row$index1]
  matrixcol2 <- count_matrix[,row$index2]
  replicate_diff <- abs(matrixcol1-matrixcol2)
  master_diff_matrix[,i] <-  replicate_diff
}


#---------------------------------------------------------------------------
# Averaging Replicates
#---------------------------------------------------------------------------
# again I'm fkinda floundering here. I bet there's better ways 

#----doing for all count tables

master_avg_matrix <- matrix(data = 0 , 
                            ncol = nrow(meta_data_group),
                            nrow = nrow(count_matrix), 
                            dimnames =  list(rownames(count_matrix), 
                                             meta_data_group$dev_stage))

for (i in 1:nrow(meta_data_group)) {
  row <- meta_data_group[i,]
  matrixcol1 <- count_matrix[,row$index1]
  matrixcol2 <- count_matrix[,row$index2]
  matrix_combined <- cbind(matrixcol1, matrixcol2)
  replicate_avg <-  round(rowMeans(matrix_combined), digits = 3)
  master_avg_matrix[,i] <-  replicate_avg
}
