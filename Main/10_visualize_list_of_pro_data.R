source(file = "Main/functions.R")
source(file = "Main/libraries.R")
source(file = "Main/global_vars.R")
source(file = "Main/03_open_expression_tables.R")
source(file = "Main/04_remove_gene_id_version.R")
source(file = "Main/05_merge.R")
source(file = "Main/09_process_l_pc.R")


#---Generate a plot to show the breakdown of our filtered expression table

filter_plots <- lapply(summary_stats, generate_filter_plot)

filter_plots
