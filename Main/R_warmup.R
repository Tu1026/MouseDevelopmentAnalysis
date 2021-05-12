## Tasks to get used to working in R and tidyverse

library(tidyverse)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)


# 1: How many unique protein coding genes are there?

glimpse(pc)
n_distinct <- pc %>%
  filter(!is.na($Symbol)%>%
  distinct(Symbol)%>%
  nrow()
print(n_distinct)


# 2: Summarize the relationship between genes and counts of transcripts (min,
# max, mean # trancripts/symbol...). Which gene has the most transcripts?

num_trans<- pc%>%
  group_by(Symbol)%>%
  summarize(num_t=n())
glimpse(pc)

pc%>%
  filter(Symbol == "")

num_trans

sum_pc <- num_trans%>%
  summarize(mean = mean(num_t), min = min(num_t), max = max(num_t))

num_trans %>%
  arrange(desc(by_group = num_t))

# 3: Describe the distribution of protein coding gene sizes. What is the average
# gene length? What is the smallest gene? Largest?