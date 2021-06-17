#This script opens all of the count tables in Data/count_tables as a list.

#Returns:
#It then saves this list as an rds object named l_expr_tables. 


source('Main/functions.R')

library(tidyverse)
library(assertthat)

#---directory where the count tables are at
count_dir <- "Data/Count_tables"


if (!(file.exists("Data/l_expr_tables.rds"))) {
  message ("l_expr_tables has not been created yet. Creating ... ")
  l_expr_tables <- open_expr_tables(count_dir)
  saveRDS(l_expr_tables, file = "Data/l_expr_tables.rds")
} else {
  message ("l_expr_tables has already been created and is sittign in /Data/")
}
