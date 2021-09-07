library(tidyverse)
library(limma)
library(edgeR)
# install.packages("gridExtra")
library(gridExtra)
contr_tests <- readRDS(file = "Data/DEAs/contr_tests.rds")

#Create table summarizing how many genes were upregulated, downregulated, and kept same, in each dev stage

count_changes <- function(column) { 

  upreg <- 0
  downreg <- 0
  same <- 0

  
  for (int in column) {
    if (int == 1) {
      upreg <- upreg + 1
    }
    if (int == -1) {
      downreg <- downreg + 1
    }
    if (int == 0) {
      same <- same + 1
    }
  }
  return (c(upreg, downreg, same))
}

counts <- data.frame(row.names = c("Up", "Down", "Same"),
                     E10.5 = count_changes(contr_tests[,1]),
                     E11.5 = count_changes(contr_tests[,2]),
                     E12.5 = count_changes(contr_tests[,3]),
                     E13.5 = count_changes(contr_tests[,4]),
                     E14.5 = count_changes(contr_tests[,5]),
                     E15.5 = count_changes(contr_tests[,6]),
                     E16.5 = count_changes(contr_tests[,7]),
                     P0 = count_changes(contr_tests[,8]))

t_counts <- as.data.frame(t(counts))

my_values <- c()
for (col in counts) { 
  my_values <- append(my_values, col)}
my_devs <- c()
for (colname in colnames(counts)) {
  my_devs <- append(my_devs, rep(colname, 3))
}

my_names <-  rep(c("up","down",'none'), 8)


df_for_bar <- data.frame(values = my_values,
                         devs = my_devs,
                         DEA_Direction = my_names)


pl_dea <- ggplot (data = df_for_bar)+
  geom_col(mapping = aes(x = devs,
                         y = values,
                         fill = DEA_Direction),
           position = "stack") + 
  labs(title = "Number of DEA Genes", x = "Dev Stage", y = "Count")
pl_dea


#--------------------------
#Sum Stats
#------------------------

updown <- append(t_counts$Up, t_counts$Down)
mean(updown)

#stage with most DEA

t_counts <- t_counts %>% 
  mutate(updown = Up+Down)
t_counts

