library(tidyverse)


top <- readRDS(file = "Data/enrichments/p_value.rds")

E10 <- head(top$allE10, 50)
E11 <- head(top$allE11, 50)
E12 <- head(top$allE12, 50)
E13 <- head(top$allE13, 50)
E14 <- head(top$allE14, 50)
E15 <- head(top$allE15, 50)
E16 <- head(top$allE16, 50)
P0 <- head(top$allP0, 50)
