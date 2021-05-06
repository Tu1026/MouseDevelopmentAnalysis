## Tasks to get used to working in R and tidyverse

library(tidyverse)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)


# 1: How many unique protein coding genes are there?

# 2: Summarize the relationship between genes and counts of transcripts (min,
# max, mean # trancripts/symbol...). Which gene has the most transcripts?

# 3: Describe the distribution of protein coding gene sizes. What is the average
# gene length? What is the smallest gene? Largest?