source(file = "Main/functions.R")
source(file = "Main/libraries.R")
source(file = "Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")

# Merge all tables in l_expr_sum with the protein coding table (pc)

# Changing the pc Gene_ID name so create_pc_expr_tables works
pc_table <- pc %>%
  rename("gene_id" = "Gene_ID")

#l_pc_expr is a list containing data frames that are our expression tables merged onto our pc
l_pc_expr <- lapply(l_expr_tables, merge_tables, pc_table=pc_table)

