library(tidyverse)
library(topGO)

# remotes::install_github("mw201608/msigdb")
library(msigdbi)

E11E10 <- readRDS(file = "Data/DEAs/E11E10.rds")

mgsa <- read.gaf(filename = "Data/mgi.gaf")
mgsa <- mgsa %>%
  dplyr::select(Symbol, GO.ID)

write_delim(x = mgsa, file = "Data/mgi.tsv", delim = "\t")

names <- unique(mgsa$Symbol)
gene2GO <- list(rep.int(0 ,times = length(names)))
myls <- vector("list", length = length(names))
names(myls) <- names

for (index in 1:nrow(mgsa)) {
  gene_symbol <- mgsa [index, "Symbol"]
  if (gene_symbol %in% names(myls)){
    myls[[gene_symbol]] <- append(myls[[gene_symbol]], mgsa[index, "GO.ID"])
  }
}

stopifnot(FALSE)

# Only get differentially expressed genes of log FC greater than abs(1)
sig <-  E11E10 %>%
  filter(adj.P.Val <= 0.05 & (logFC >=1 | logFC <=-1))
pos_sig <- E11E10 %>%
  filter(adj.P.Val <= 0.05 & logFC >=1)
neg_sig <- E11E10 %>%
  filter(adj.P.Val <= 0.05 & logFC <=-1)

# Coerce all diff expressed genes into list of gene_name:logFC
list_sig <- as.list(sig$logFC)
names(list_sig) <- rownames(sig)

#actually. I think we shouldn't prematurely only use the sig genes. So i'll coerce the origonal as well
vector_all <- as_vector(E11E10$adj.P.Val)
names(vector_all) <- rownames(E11E10)

#only cerate selection criteria

selection_mechanism <- function(list_sig){
  #select ALL genes. This assumes you've pre-selected your genes of interest
  output_list <- list()
  for (gene in names(list_sig)) {
    output_list[gene] = TRUE
  }

  return(output_list)
}


#Creating godata for the annot function
my_gene2GO <- as.vector(mgsa$GO.ID)
names(my_gene2GO) <- mgsa$Symbol


# Create GOdata object
GOdata <- new("topGOdata",
              description = "E11E10 topGO",
              ontology = "BP",
              allGenes = vector_all,
              geneSelectionFun = topDiffGenes,
              nodeSize = 5,
              annot = annFUN.gene2GO, gene2GO = my_gene2GO)
              


#Accessing GoData
description(GOdata)

genes <- genes(GOdata)
numGenes(GOdata)

selGenes <- sample(genes, 10)
selGenes                         

gs <- geneScore(GOdata, whichGenes = selGenes)
gs

vector_all[selGenes]
topDiffGenes(selGenes)

sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)


graph(GOdata)

ug <- usedGO(GOdata)
ug


sel.terms <- sample(usedGO(GOdata),10)
num_ann_genes <- countGenesInTerm(GOdata, sel.terms)
num_ann_genes
ann.genes <- genesInTerm(GOdata, sel.terms)
ann.genes


ann.score <- scoresInTerm(GOdata, sel.terms)
ann.score

term.stats <-  termStat(GOdata, sel.terms)
term.stats


# Running the enrichment tests

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFis

pvalFis <- score(resultFis)
head(pvalFis)


hist(pvalFis, 50, xlab = "p-values")

geneData(resultFis)


allRes <- GenTable(GOdata,
                   classic = resultFis,
                   orderBy = classic,
                   topNodes = length(score(resultFis)))
                   

goID <- allRes[99, "GO.ID"]
print(showGroupDensity(GOdata, goID, ranks = TRUE))

BiocManager::install("Rgraphviz")
library(Rgraphviz)

a <- showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
a$complete.dag
a$dag
