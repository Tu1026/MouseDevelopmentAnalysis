library(tidyverse)
library(topGO)

# remotes::install_github("mw201608/msigdb")
library(msigdbi)

# BiocManager::install("Rgraphviz")
library(Rgraphviz)


# Load the GO ontology mapping
if (! "gene2GO" %in% ls()) {
  gene2GO <- readRDS(file = "Data/GO_mappings.rds")
  }


# test1 <- ontology_map %>%
#   filter(GO.ID == "GO:0030900")
# 
# built <- c()
# 
# for (index in 1:length(my_gene2GO)){
#   if (my_gene2GO[[index]] == "GO:0030900") {
#     built <- append(my_gene2GO[index], built)
#   }
# }
# 
# test1$Symbol %in% names(built)
# 
# built2 <- c()
# 
# for (index in 1:length(my_gene2GO)){
#   if (names(my_gene2GO[index]) == "Sox1") {
#     built2 <- append(my_gene2GO[index], built2)
#   }
# }

desired_DEAs <- c("Data/DEAs/allE10.rds",
                  "Data/DEAs/allE11.rds",
                  "Data/DEAs/allE12.rds",
                  "Data/DEAs/allE13.rds",
                  "Data/DEAs/allE14.rds",
                  "Data/DEAs/allE15.rds",
                  "Data/DEAs/allE16.rds",
                  "Data/DEAs/allP0.rds")
          

# Create selection criteria for logFC GO enrichment analysis
logFC_upregulate <- function(allScore){
  #select ALL genes. This assumes you've pre-selected your genes of interest
  return (allScore >= 1)
}
logFC_downregulate <- function(allScore){
  #select ALL genes. This assumes you've pre-selected your genes of interest
  return (allScore <= -1)
}
low_p <- function(allScore) {
  return( allScore < 0.01 )
}


do_enrichment <- function(DEA.rds, my_gene2GO, ontology_type, node_size, selection_method) {
  # given a DEA, use topGO to perform an GO enrichment analysis for ONLY the differentially expressed genes
  #with respect to logFC
  
  #@PARAM:
  #my_gene2GO: A vector containing genename:GO mappings
  #ontology_type: The type of GO ontology you want to use. Can be BP, MF, or CC
  #node_size: The min amount of gene to ontology mappings that are required for the mapping to be considered significant
  #selection_method: How the enrichment selects for GO terms that are enriched. logFC_upregulate will enrich for upregulated genes, whereas
  #  logFC_downregulate will enrich for downregulated genes

  #open target DEA rds object
  DEA <- readRDS(file = DEA.rds)
  
  #quickly extract the name of the DEA just for naming purposes later on
  DEA_name <- str_replace(string = str_split(string = DEA.rds, pattern = "/")[[1]][3], pattern = ".rds", replacement = "")
  message(DEA_name)

  ####- The commented out section was used when doing logFC
  # # Only get differentially exprepressed genes with alpha 0.01
  # sig <-  DEA %>%
  #   filter(adj.P.Val <= 0.05)
  # if (is.null(sig)) {
  #   warning(paste(DEA_name, "did not have any DEA"))
  #   return(NA)
  # }
  # 
  # # Coerce all diff expressed genes into vector of gene_name:logFC
  # my_allGenes <- as.vector(sig$logFC)
  # names(my_allGenes) <- rownames(sig)

  my_allGenes <- as.vector(DEA$adj.P.Val)
  names(my_allGenes) <- rownames(DEA)
  
  
  
  # Create GOdata object
  GOdata <- new("topGOdata",
                description = paste0(DEA_name,"_p"),
                ontology = ontology_type,
                allGenes = my_allGenes,
                geneSelectionFun = selection_method,
                nodeSize = node_size,
                annot = annFUN.gene2GO,
                gene2GO = my_gene2GO) 
  
  # Run fisher test
  
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  

  if (length(score(resultFis)) == 0){
    warning(paste(DEA_name, "did not have any significant annotations"))
    return(NA)
  }
  
  # Create ResultTop object

  
  allRes <- GenTable(GOdata,
                     classic = resultFis,
                     orderBy = classic,
                     topNodes = length(score(resultFis)))

  return(allRes)
}


DEA_names <- lapply(desired_DEAs, function(x) str_replace(string = str_split(string = x, pattern = "/")[[1]][3], pattern = ".rds", replacement = ""))

p_value_enr <- lapply(desired_DEAs, do_enrichment,
               my_gene2GO = gene2GO,
               ontology_type = "BP",
               node_size = 10,
               selection_method =low_p)

names(p_value_enr) <- DEA_names

saveRDS(p_value_enr, file = "Data/enrichments/p_value.rds")

stopifnot(FALSE)


BPenrichments_downregulate <- lapply(desired_DEAs, do_enrichment,
                                 my_gene2GO = gene2GO,
                                 ontology_type = "BP",
                                 node_size = 10,
                                 selection_method =logFC_downregulate)
names(BPenrichments_downregulate) <- DEA_names
# 
# MFenrichments_upregulate <- lapply(desired_DEAs, do_enrichment,
#                                    my_gene2GO = my_gene2GO,
#                                    ontology_type = "MF",
#                                    node_size = 10,
#                                    selection_method =logFC_upregulate)
# names(MFenrichments_upregulate) <- DEA_names
# 
# CCenrichments_upregulate <- lapply(desired_DEAs, do_enrichment,
#                                    my_gene2GO = my_gene2GO,
#                                    ontology_type = "CC",
#                                    node_size = 10,
#                                    selection_method =logFC_upregulate)
# names(CCenrichments_upregulate) <- DEA_names

dir.create("Data/enrichments", showWarnings = FALSE)
saveRDS(BPenrichments_upregulate, file = "Data/enrichments/BPenrichments_upregulate.rds")
saveRDS(BPenrichments_downregulate, file = "Data/enrichments/BPenrichments_downregulate.rds")


#TODO. Create all possible upregu and downreg for the 3 GO types. Save them
#Then start to look at specific GO's like GO:0030900 - forbrain development







# # Only get differentially expressed genes of log FC greater than abs(1)
# sig <-  DEA %>%
#   filter(adj.P.Val <= 0.01)
# 
# # Coerce all diff expressed genes into list of gene_name:logFC
# my_allGenes <- as.vector(sig$logFC)
# names(my_allGenes) <- rownames(sig)
# 
# # Create selection criteria for locFC
# 
# logFC_select <- function(allScore){
#   #select ALL genes. This assumes you've pre-selected your genes of interest
#   return (allScore >= 1 | allScore <= -1)
# }
# 
# #Creating godata for the annot function
# # so it can link the gene symbols to Go annotations
# my_gene2GO <- as.vector(mgsa$GO.ID)
# names(my_gene2GO) <- mgsa$Symbol
# 
# 
# # Create GOdata object
# GOdata <- new("topGOdata",
#               description = "E11E10 topGOFC",
#               ontology = "BP",
#               allGenes = my_allGenes,
#               geneSelectionFun = logFC_select,
#               nodeSize = 10,
#               annot = annFUN.gene2GO, gene2GO = my_gene2GO)
# 
# 
# 
# #Accessing GoData
# description(GOdata)
# 
# genes <- genes(GOdata)
# numGenes(GOdata)
# 
# selGenes <- sample(genes, 10)
# selGenes                         
# 
# gs <- geneScore(GOdata, whichGenes = selGenes)
# gs
# 
# my_allGenes[selGenes]
# logFC_select(selGenes)
# 
# sg <- sigGenes(GOdata)
# str(sg)
# numSigGenes(GOdata)
# 
# 
# graph(GOdata)
# 
# ug <- usedGO(GOdata)
# ug
# 
# 
# sel.terms <- sample(usedGO(GOdata),10)
# num_ann_genes <- countGenesInTerm(GOdata, sel.terms)
# num_ann_genes
# ann.genes <- genesInTerm(GOdata, sel.terms)
# ann.genes
# 
# 
# ann.score <- scoresInTerm(GOdata, sel.terms)
# ann.score
# 
# term.stats <-  termStat(GOdata, sel.terms)
# term.stats
# 
# 
# # Running the enrichment tests
# 
# resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
# resultFis
# 
# pvalFis <- score(resultFis)
# head(pvalFis)
# 
# 
# hist(pvalFis, 50, xlab = "p-values")
# 
# geneData(resultFis)
# 
# 
# allRes <- GenTable(GOdata,
#                    classic = resultFis,
#                    orderBy = classic,
#                    topNodes = length(score(resultFis)))
# 
# view(allRes)
# goID <- allRes[99, "GO.ID"]
# print(showGroupDensity(GOdata, goID, ranks = TRUE))
# 
# 
# 
# a <- showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
# a$complete.dag
# a$dag
