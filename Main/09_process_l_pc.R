#
# Do 06, 07, 08. But do it with l_pc_expr using lapply
#

source(file = "Main/functions.R")
source(file = "Main/libraries.R")
source(file = "Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")
source(file = "Main/05_merge.R")

summary_stats <- lapply(l_pc_expr, generate_summary_stats)
  

