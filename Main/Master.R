#
# This script will retnrun you the l_Expressed_Genes, and its affiliated metadata: replicate_meta.
#
source(file = "Main/functions.R")
source(file = "Main/libraries.R")
source(file = "Main/global_vars.R")
source('Main/02_create_metadata.R')
source('Main/11_select_largest_transcripts_from_list.R')

print("yay")

glimpse(l_Expressed_Genes)
glimpse(replicate_meta)
