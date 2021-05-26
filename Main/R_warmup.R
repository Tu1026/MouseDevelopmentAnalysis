## Tasks to get used to working in R and tidyverse

library(tidyverse)

pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)


# 1: How many unique protein coding genes are there?

glimpse(pc)
n_distinct <- pc %>%
  filter(pc$Symbol != "")%>%
  distinct(Symbol)%>%
  nrow()
print(n_distinct)
21786

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

num_trans <-  num_trans %>%
  filter(num_trans$Symbol !="")%>%
  arrange(desc(by_group = num_t))

print(paste(num_trans$Symbol[1],"has the most transcripts"))

# 3: Describe the distribution of protein coding gene sizes. What is the average
# gene length? What is the smallest gene? Largest?

glimpse(pc)
lengths <- pc %>%
  select(Symbol,Start,End)%>%
  mutate(length = End-Start)%>%
  arrange(desc(by_group = length))

largest_gene <- lengths[1,]
smallest_gene <- lengths[nrow(lengths),]

length_summary <- lengths%>%
  summarize(mean = mean(lengths$length), min = min(lengths$length), max = max(lengths$length))


