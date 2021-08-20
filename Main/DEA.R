library(tidyverse)
library(limma)
library(edgeR)

meta <- read.delim("Data/complete_meta_data.tsv", stringsAsFactors = FALSE)
pc <- read.delim("Data/ensembl_mouse_protein_coding_104.tsv", stringsAsFactors = FALSE)
expected_mat <- readRDS(file = "Data/pc_count_matrix_expected.rds")

#--------------------------------------------------------------
# Normalization and Filtering
#--------------------------------------------------------------



#Create a DGEList object using edgeR: DGE is a digital gene expression data class
#its just an S4 class for storing read counts

dge <- DGEList(counts = expected_mat)



#Visualizing removing genes
#convert our  values to cpm and logCPM values so that we can get some visuzliation
cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE, prior.count=3)

#get mean and median library size counts for visualization
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6

#Visualize adding non-filtered data to plot
lcpm.cutoff <- log2(10/M + 2/L)
samplenames <- paste0(colnames(dge),"_",meta$dev_stage)
library(RColorBrewer)
nsamples <- ncol(dge)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

#Remove low counts
#filterByExpr is a EdgeR fuinction that keeps genes that have at least a certain min value
# in a sufficient number of samples. In other words, the filtering keeps genes that have count-per-million above
# k in n samples. where k is determined by mincount and the librazy sizes and n is determined by the design matrix.
# filterByExpr automatically converts our dge to CPM values.
groups <- factor(meta$dev_stage)
keep <- filterByExpr(dge, groups)
dge <- dge[keep, keep.lib.sizes = FALSE]


# -- ADD filtered data to plot
# legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(dge, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
# legend("topright", samplenames, text.col=col, bty="n")


#See some removed
removed <- expected_mat[!keep,]
removedsum <- rowSums(removed)
removedavg <- rowMeans(removed)


#to TMM Normalization (trimmed mean of M-values normalization method)
#normalizes for batch effects(depth) and gene length
dge <- calcNormFactors(dge)


# Unsupervised clustering of samples -already done with PCA. should be fine


#--------------------------------------------------------------
# Limma DEA
#--------------------------------------------------------------

#Create Design
design <- model.matrix( ~ 0+groups)
colnames(design) <- levels(groups)
rownames(design) <- meta$id

#Setting up Contrasts
# we need to tell limma what groups we want to look for DEA in - sett up the pairwise comparisons

allcontrasts <- list(  first3 = (E11.5+E10.5+E12.5)/3 - (E13.5 + E14.5 + E15.5+ E16.5 + P0)/5,
                       allP0 =  P0    - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5)/7,
                       allE10 = E10.5 - (E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
                       allE11 = E11.5 - (E10.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
                       allE12 = E12.5 - (E10.5 + E11.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
                       allE13 = E13.5 - (E10.5 + E11.5 + E12.5 + E14.5 + E15.5 + E16.5 + P0)/7,
                       allE14 = E14.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E15.5 + E16.5 + P0)/7,
                       allE15 = E15.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E16.5 + P0)/7,
                       allE16 = E16.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + P0)/7,
                       
                       E11E10 = E11.5 - E10.5,
                       E12E11 = E12.5 - E11.5,
                       E13E12 = E13.5 - E12.5,
                       E14E13 = E14.5 - E13.5,
                       E15E14 = E15.5 - E14.5,
                       E16E15 = E16.5 - E15.5,
                       P0E16  = P0    - E16.5)

contrast <- makeContrasts(
  first3 = (E11.5+E10.5+E12.5)/3 - (E13.5 + E14.5 + E15.5+ E16.5 + P0)/5,
  allP0 =  P0    - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5)/7,
  allE10 = E10.5 - (E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
  allE11 = E11.5 - (E10.5 + E12.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
  allE12 = E12.5 - (E10.5 + E11.5 + E13.5 + E14.5 + E15.5 + E16.5 + P0)/7,
  allE13 = E13.5 - (E10.5 + E11.5 + E12.5 + E14.5 + E15.5 + E16.5 + P0)/7,
  allE14 = E14.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E15.5 + E16.5 + P0)/7,
  allE15 = E15.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E16.5 + P0)/7,
  allE16 = E16.5 - (E10.5 + E11.5 + E12.5 + E13.5 + E14.5 + E15.5 + P0)/7,

  E11E10 = E11.5 - E10.5,
  E12E11 = E12.5 - E11.5,
  E13E12 = E13.5 - E12.5,
  E14E13 = E14.5 - E13.5,
  E15E14 = E15.5 - E14.5,
  E16E15 = E16.5 - E15.5,
  P0E16  = P0    - E16.5,

  levels = colnames(design)
)


# Removing heteroscedascity from count data
# Methods that model the heteroscedasity use a negative binomial distribution assume a quadratic mean-variance relationship
# In limma, we simply linear model using log-CPM values- which are heteroscetdatic,
# but to acount for it, the mean-variance relationships is accomodated using precision weights calculated by the voom function

# mea_variance plot can be easily calculated
v <- voom(dge, design, plot=TRUE)
v
#The fact that the line increase at the beggining suggests that small reads,which usually have low,
# were mostly removed. This suggets that our level of filtering was sufficient
# we will have to add our fitted model later onto this using


# Fitting linear models for comparisons of interest
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contrast)
efit <- eBayes(vfit)
plotSA(efit, main = "Final model: Mean-variance trend")


# Examining the number of DE genes


# Examining individual DE genes from top to bottom
#use topTable()
topTable(efit)

# saveRDS(object = efit, file = "Data/DEA.rds")

