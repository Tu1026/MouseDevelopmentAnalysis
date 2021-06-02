
library(tidyverse)

libraries <- c("assertthat", "biomaRt", "tools", "googlesheets4", "tools")

for (i in libraries) {
  if (!(i %in% library(quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)$results[,1])) {
    library(i)
  }
}

