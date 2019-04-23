####################################################################################
# Installing Packages in R and Running an RNA Sequence Analysis
###################################################################################
#
# Resources for understanding RNA Seq
#
# RNA-sqlopedia:  https://rnaseq.uoregon.edu/ 
#
# Workflow Example: 
#
#  https://www.bioconductor.org/help/workflows/RNAseq123/ 
#  https://master.bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# 
#  Step 0: 
#          from the command line install libmariadb-client-lgpl-dev:
#          sudo apt-get update
#          sudo apt-get install libmariadb-client-lgpl-dev
#          sudo apt-get install libcurl4-openssl-dev libssl-dev
#
#          You must run the installation portion of this script as an administrative
#          user because it compiles and installs packages!
#  
#  Step 1: If not already done correctly, re-start rstudio as administrator (we want to install packages.)
#          from the command line type: sudo rstudio
#
#
# Get R utils and tell R we want packages from bionconductor
#

install.packages("R.utils")
source("https://bioconductor.org/biocLite.R")

# Install Limma, Glimma, edgeR and us.musculus Packages:
# run each of these one at a time - only need to install the very first time
# Rerun these if you update RStudio or Bioconductor

#biocLite("limma")
#biocLite("Glimma")
#biocLite("edgeR")
#biocLite("openssl")
#biocLite("Mus.musculus")

# on subsequent runs 

library(limma)
library(Glimma)
library(edgeR)
# Transcriptome profiling of purified mouse mammary stem, progenitor and mature cell populations
library(Mus.musculus)  

#
# Step 2: Make a folder to store files for the project
#

dir.create("RNA_Seq")  # Create a folder on the desktop to hold the RNA_Seq data
setwd("./RNA_Seq")     # change current folder in R to the new folder you just created

#
# Step 3: Read in the counts data for this RNA Seq tutorial
#
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"  # set the URL we will get stuff from
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")            # download the tar file
utils::untar("GSE63310_RAW.tar", exdir = ".")                                # extract all the files
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

# Each of these text files contains the raw gene-level counts for a given sample.
# Note that our analysis only includes the basal, LP and ML samples from this experiment 
# (see associated file names below).

files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", 
           "GSM1545538_purep53.txt", "GSM1545539_JMS8-2.txt", 
           "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", 
           "GSM1545545_JMS9-P8c.txt")
read.delim(files[1], nrow=5)

# Whilst each of the nine text files can be read into R separately and combined 
# into a matrix of counts, edgeR offers a convenient way to do this in one step 
# using the readDGE function. The resulting DGEList-object contains a matrix of 
# counts with 27,179 rows associated with unique Entrez gene identifiers (IDs) and 
# nine columns associated with the individual samples in the experiment.

# Store all the data from the columns into our data frame called x
x <- readDGE(files, columns=c(1,3))
class(x)
dim(x)

#####################################################################################################
# 4.2 Organising sample information
#####################################################################################################

# For downstream analysis, sample-level information related to the experimental design needs
# to be associated with the columns of the counts matrix. This should include experimental variables, 
# both biological and technical, that could have an effect on expression levels. 
#
# The cell types for this experiment:  
#
#     basal - a type of cell in the innermost layer of the epidermis or other epithelial tissue.
#     LP - Luminal progenitor-enriched and 
#     ML - mature luminal-enriched 
#     see following link for more details: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310
#
# The adult mammary epithelium is composed of basal and luminal cells. The luminal lineage comprises two major 
# cell populations, positive and negative for estrogen and progesterone receptors (ER and PR, respectively), 
# both containing clonogenic progenitor cells (all cells are clones of one another).
#
# Other variables are:
#    genotype (wild-type, knock-out), 
#    phenotype (disease status, sex, age), 
#    sample treatment (drug, control) and 
#    batch information (date experiment was performed if samples were collected and analysed at distinct 
#                       time points) to name just a few.
# Our DGEList-object contains a samples data frame that stores both cell type (or group) 
# and batch (sequencing lane) information, each of which consists of three distinct levels. 
# Note that within x$samples, library sizes are automatically calculated for each sample and 
# normalisation factors are set to 1. For simplicity, we remove the GEO sample IDs (GSM*) 
# from the column names of our DGEList-object x.

samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames

# set all the columnames to our sample names in our data frame x
colnames(x) <- samplenames

# create a factor array of cell types for each column
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))

# store the groups in our data frame x as $samples$groups
x$samples$group <- group

# when rna sequencing was run, it was run on multiple lanes. Create a factor list that 
# identifies the lanes and store this in our x data set as x$samples$lane
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane

#observe the samples and labeling information for each sample
x$samples

###################################################################################
# 4.3 Organizing Gene Annotations
###################################################################################

# A second data frame named genes in the DGEList-object is used to store gene-level 
# information associated with rows of the counts matrix. 
# This information can be retrieved using organism specific packages such as Mus.musculus 
# (Bioconductor Core Team 2016b) for mouse 
# (or Homo.sapiens (Bioconductor Core Team 2016a) for human) or the 
# (biomaRt package Durinck et al. 2005, 2009) which interfaces the Ensembl genome databases 
# in order to perform gene annotation. 
#
# About gene annotations
#
# The type of gene annotation information that can be retrieved includes:
#   gene symbols, 
#   gene names, 
#   chromosome names and locations, 
#   Entrez gene IDs - (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3013746/ ) 
#   Refseq gene IDs  - (see https://www.ncbi.nlm.nih.gov/refseq/ )
#   Ensembl gene IDs - ( see http://www.ensembl.org/Multi/Search/Results?q=ENSG00000010404%3Fredirect%3Dno )
#   are just a few of the databases of gene annotations out there... 
#   biomaRt primarily works off Ensembl gene IDs, 
#   Mus.musculus packages information from various sources and allows users to choose
#   between many different gene IDs as the key. 
#   The Entrez gene IDs available in our dataset were annotated using the Mus.musculus package 
#   to retrieve associated gene symbols and chromosome information.

# Set a vector to the geneid's of x
geneid <- rownames(x)

# Look at just the first few geneids
head(geneid)

# Now create a mapping structure that will map gene IDs (if possible) to Entrez gene IDs 
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")

# Look at the first few entries of our geneid mapping tables
head(genes)

# As with any gene ID, Entrez gene IDs may not map one-to-one to the gene information of interest. 
# It is important to check for duplicated gene IDs and to understand the source of duplication before
# resolving them. 
#
# Our gene annotation contains 28 genes that map to multiple chromosomes (e.g. gene Gm1987 
# is associated with chr4 and chr4_JH584294_random and microRNA Mir5098 is associated with 
# chr2, chr5, chr8, chr11 and chr17). 
#
# To resolve duplicate gene IDs one could combine all chromosome information from the multi-mapped genes, 
# such that gene Gm1987 would be is assigned to chr4 and chr4_JH584294_random, or select one of the 
# chromosomes to represent the gene with duplicate annotation. 
#
# For simplicity we do the latter, keeping only the first occurrence of each gene ID.
#
genes <- genes[!duplicated(genes$ENTREZID),]

# In this example, the gene order is the same in both the annotation and the data object. 
# If this is not the case due to missing and/or rearranged gene IDs, the match function can be used to 
# order genes correctly. 
#
# The data frame of gene annotations is then added to the data object and neatly packaged in a 
# DGEList-object containing raw count data with associated sample information and gene annotations.
#

x$genes <- genes
x

###################################################################################
# 5 Data Pre-processing
###################################################################################

###################################################################################
# 5.1 Transformations from the raw-scale
###################################################################################

# For differential expression and related analyses, gene expression is rarely considered at the 
# level of raw counts since libraries sequenced at a greater depth will result in higher counts. 
# Rather, it is common practice to transform raw counts onto a scale that accounts for such library 
# size differences. 
#
# Popular transformations include:
# counts per million (CPM), 
# log2-counts per million (log-CPM), 
# reads per kilobase of transcript per million (RPKM), and 
# fragments per kilobase of transcript per million (FPKM).
#
# In our analyses, CPM and log-CPM transformations are used regularly although they do not account 
# for gene length differences as RPKM and FPKM values do. Whilst RPKM and FPKM values can just as well 
# be used, CPM and log-CPM values can be calculated using a counts matrix alone and will suffice for 
# the type of comparisons we are interested in. 
# 
# Assuming that there are no differences in isoform usage between conditions, 
# differential expression analyses look at gene expression changes between conditions 
# rather than comparing expression across multiple genes or drawing conclusions on absolute 
# levels of expression. In other words, gene lengths remain constant for comparisons of interest 
# and any observed differences are a result of changes in condition rather than changes in gene length.
#
# Here raw counts are converted to CPM and log-CPM values using the cpm function in edgeR. 
# RPKM values are just as easily calculated as CPM values using the rpkm function in edgeR 
# if gene lengths are available.

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# A CPM value of 1 for a gene equates to having 20 counts in the sample with the lowest sequencing
# depth (JMS0-P8c, library size approx. 20 million) or 76 counts in the sample with the greatest 
# sequencing depth (JMS8-3, library size approx. 76 million).
# 
# The log-CPM values will be used for exploratory plots. When log=TRUE, the cpm function adds an 
# offset to the CPM values before converting to the log2-scale. 
# 
# By default, the offset is 2/L where 2 is the “prior count” and L is the average library size in 
# millions, so the log-CPM values are related to the CPM values by log2(CPM + 2/L). 
# This calculation ensures that any two read counts with identical CPM values will also have identical 
# log-CPM values. 
# 
# The prior count avoids taking the logarithm of zero, and also reduces spurious variability for 
# genes with very low counts by shrinking all the inter-sample log-fold-changes towards zero, 
# something that is helpful for exploratory plotting. 
# 
# For this dataset, the average library size is about 45.5 million, so L approx. 45.5 and the 
# minimum log-CPM value for each sample becomes log2(2/45.5) = -4.51. In other words, a count 
# of zero for this data maps to a log-CPM value of -4.51 after adding the prior count or offset:

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

# Log-CPM values are also used in downstream linear modeling via limma’s voom function, although 
# voom recomputes its own log-CPM values internally with a smaller prior count.

###################################################################################
# 5.2 Removing genes that are lowly expressed
###################################################################################

# All datasets will include a mix of genes that are expressed and those that are not expressed. 
# Whilst it is of interest to examine genes that are expressed in one condition but not in another, 
# some genes are unexpressed throughout all samples. In fact, 19% of genes in this dataset have zero 
# counts across all nine samples.

table(rowSums(x$counts==0)==9)

# Plotting the distribution log-CPM values shows that a sizeable proportion of genes within 
# each sample are either unexpressed or lowly-expressed with log-CPM values that are 
# small or negative.
#
# Genes that do not have a worthwhile number of reads in any sample should be filtered out of 
# the downstream analyses. There are several reasons for this. From a biological point of view, 
# genes that not expressed at a biologically meaningful level in any condition are not of interest 
# and are therefore best ignored. From a statistical point of view, removing low count genes allows 
# he mean-variance relationship in the data to be estimated with greater reliability and also reduces 
# the number of statistical tests that need to be carried out in downstream analyses looking at differential 
# expression.
# 
# The filterByExpr function in the edgeR package provides an automatic way to filter genes, while keeping 
# as many genes as possible with worthwhile counts.
#
# By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, 
# where the number of samples is chosen according to the minimum group sample size. The actual filtering 
# uses CPM values rather than counts in order to avoid giving preference to samples with large library sizes. 
# For this dataset, the median library size is about 51 million and 10/51 approx. 0.2, so the filterByExpr 
# function keeps genes that have a CPM of 0.2 or more in at least three samples. A biologically interesting 
# gene should be expressed in at least three samples because all the cell type groups have three replicates. 
# The cutoffs used depend on the sequencing depth and on the experimental design. If the library sizes had been 
# larger then a lower CPM cutoff would have been chosen, because larger library sizes provide better resolution 
# to explore more genes at lower expression levels. 
#
# Alternatively, smaller library sizes decrease our ability to explore marginal genes and hence would have led to 
# a higher CPM cutoff.

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# Using this criterion, the number of genes is reduced to 16,624, about 60% of the number that we started with 
# (panel B of the next figure). Note that subsetting the entire DGEList-object removes both the counts and the 
# associated gene information for the filtered genes. The filtered DGEList-object keeps the gene information and
# the counts for the retained genes correctly associated.
#
# Code to produce the figure is given below.
#
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#
# Figure 1: The density of log-CPM values for raw pre-filtered data (A) 
# and post-filtered data (B) are shown for each sample
# Dotted vertical lines mark the log-CPM threshold (equivalent to a CPM value 
# of about 0.2) used in tBasal_Lhe filtering step.
#

###################################################################################
# 5.3 Normalising gene expression distributions
###################################################################################

# During the sample preparation or sequencing process, external factors that are not 
# of biological interest can affect the expression of individual samples. 
# For example, samples processed in the first batch of an experiment can have higher 
# expression overall when compared to samples processed in a second batch. 
# It is assumed that all samples should have a similar range and distribution of expression 
# values. Normalisation is required to ensure that the expression distributions of each sample 
# are similar across the entire experiment.
#
# (AN ASSIDE NOT IN THE ORIGINAL PAPER) there are those who feel that this normalization step is not accurate, because different
# experiments may have different distibutions and chacteristics. This is an open area of research.
#
# Any plot showing the per sample expression distributions, such as a density or boxplot, is useful 
# in determining whether any samples are dissimilar to others. Distributions of log-CPM values are 
# similar throughout all samples within this dataset (panel B of the figure above).
# 
# Nonetheless, normalisation by the method of trimmed mean of M-values (TMM) (Robinson and Oshlack 2010) 
# is performed using the calcNormFactors function in edgeR. The normalisation factors calculated here are 
# used as a scaling factor for the library sizes. When working with DGEList-objects, these normalisation 
# factors are automatically stored in  x$samples$norm.factors. 
#
# For this dataset the effect of TMM-normalisation is mild, as evident in the magnitude of the scaling factors, 
# which are all relatively close to 1.
#

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

# To give a better visual representation of the effects of normalisation, the data was duplicated then 
# adjusted so that the counts of the first sample are reduced to 5% of their original values, and in the 
# second sample they are inflated to be 5-times larger.
# 

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

# The figure below shows the expression distribution of samples for unnormalised and normalised data, 
# where distributions are noticeably different pre-normalisation and are similar post-normalisation. 
# Here the first sample has a small TMM scaling factor of 0.06, whereas the second sample has a large scaling factor of 6.08
# – neither values are close to 1.
# 

par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)Basal_L
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#
# Figure 2: Example data: Boxplots of log-CPM values showing expression distributions for 
# unnormalised data (A) and normalised data (B) for each sample in the modified dataset where the counts in 
# samples 1 and 2 have been scaled to 5% and 500% of their original values respectively
#

###################################################################################
# 5.4 Unsupervised clustering of samples
###################################################################################

# In our opinion, one of the most important exploratory plots to examine for gene expression analyses is the
# multi-dimensional scaling (MDS) plot, or similar. The plot shows similarities and dissimilarities between samples 
# in an unsupervised manner so that one can have an idea of the extent to which differential expression can be detected 
# before carrying out formal tests. Ideally, samples would cluster well within the primary condition of interest, and any 
# sample straying far from its group could be identified and followed up for sources of error or extra variation. If present, 
# technical replicates should lie very close to one another.
#
# Such a plot can be made in limma using the plotMDS function. The first dimension represents the leading-fold-change that 
# best separates samples and explains the largest proportion of variation in the data, with subsequent dimensions having a 
# smaller effect and being orthogonal to the ones before it. 
#
# When experimental design involves multiple factors, it is recommended that each factor is examined over several dimensions. 
# If samples cluster by a given factor in any of these dimensions, it suggests that the factor contributes to expression 
# differences and is worth including in the linear modelling. On the other hand, factors that show little or no effect may 
# be left out of downstream analysis.
#
# In this dataset, samples can be seen to cluster well within experimental groups over dimension 1 and 2, and then separate by 
# sequencing lane (sample batch) over dimension 3 (shown in the plot below). Keeping in mind that the first dimension explains 
# the largest proportion of variation in the data, notice that the range of values over the dimensions become smaller as we move 
# to higher dimensions.
#
# Whilst all samples cluster by groups, the largest transcriptional difference is observed between basal and LP, 
# and basal and ML over dimension 1. For this reason, it is expected that pairwise comparisons between cell populations 
# will result in a greater number of DE genes for comparisons involving basal samples, and relatively small numbers of 
# DE genes when comparing ML to LP. Datasets where samples do not cluster by experimental group may show little or no evidence 
# of differential expression in the downstream analysis.
#
# To create the MDS plots, we assign different colours to the factors of interest. Dimensions 1 and 2 are examined using 
# the color grouping defined by cell types.
#
# Dimensions 3 and 4 are examined using the colour grouping defined by sequencing lanes (batch).
#
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

# Figure 3: MDS plots of log-CPM values over dimensions 1 and 2 with samples coloured 
# and labeled by sample groups (A) and over dimensions 3 and 4 with samples coloured 
# and labeled by sequencing lane (B)
# Distances on the plot correspond to the leading fold-change, which is the average 
# (root-mean-square) log2-fold-change for the 500 genes most divergent between each 
# pair of samples by default.
#

# Alternatively, the Glimma package offers the convenience of an interactive 
# MDS plot where multiple dimensions can be explored. 
# The glMDSPlot function generates an html page (that is opened in a browser 
# if launch=TRUE) with an MDS plot in the left panel and a barplot showing the 
# proportion of variation explained by each dimension in the right panel. 
# Clicking on the bars of the bar plot changes the pair of dimensions plotted 
# in the MDS plot, and hovering over the individual points reveals the sample label. 
# The colour scheme can be changed as well to highlight cell population or sequencing 
# lane (batch). An interactive MDS plot of this 
# dataset can be found at
# http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MDS-Plot.html.
#

glMDSPlot(lcpm, labels=paste(group, lane, sep="_"),
          groups=x$samples[,c(2,5)], launch=FALSE)


###################################################################################
# 6 Differential expression analyasis
###################################################################################

###################################################################################
# 6.1 Creating a design matrix and contrasts
###################################################################################

# In this study, it is of interest to see which genes are expressed at 
# different levels between the three cell populations profiled. 
#
# In our analysis, linear models are fitted to the data with the assumption 
# that the underlying data is normally distributed. 
#
# To get started, a design matrix is set up with both the cell population and sequencing 
# lane (batch) information.

# for understanding design matrices in R you might want to read:
# https://genomicsclass.github.io/book/pages/expressing_design_formula.html 
# 

# Here they are removing the intercept from group but leaving an intercept for lane.
#

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

# For a given experiment, there are usually several equivalent ways to set up an appropriate 
# design matrix. For example, ~0+group+lane removes the intercept from the first factor, group, 
# but an intercept remains in the second factor lane. 
#
# Alternatively, ~group+lane could be used to keep the intercepts in both group and lane. 
# Understanding how to interpret the coefficients estimated in a given model is key here. 
# We choose the first model for our analysis, as setting up model contrasts is more straight 
# forward in the absence of an intercept for group. 
#
# Contrasts for pairwise comparisons between cell populations are set up in limma using 
# the makeContrasts function.
#

contr.matrix <- makeContrasts(
  BasalvsLP = Basal-LP, 
  BasalvsML = Basal - ML, 
  LPvsML = LP - ML, 
  levels = colnames(design))
contr.matrix

# A key strength of limma’s linear modelling approach, is the ability accommodate arbitrary experimental complexity. 
# Simple designs, such as the one in this workflow, with cell type and batch, through to more complicated 
# factorial designs and models with interaction terms can be handled relatively easily. 
# Where experimental or technical effects can be modelled using a random effect, another possibility in limma 
# is to estimate correlations using duplicateCorrelation by specifying a block argument for both this function 
# and in the lmFit linear modelling step.
#

###################################################################################
# 6.2 Removing heteroscedascity from count data
###################################################################################

# It has been shown that for RNA-seq count data, the variance is not independent of the mean 
# (Law et al. 2014) – this is true of raw counts or when transformed to log-CPM values. Methods 
# that model counts using a Negative Binomial distribution assume a quadratic mean-variance relationship.
# See https://en.wikipedia.org/wiki/Negative_binomial_distribution 
#
# In limma, linear modelling is carried out on the log-CPM values which are assumed to be normally distributed 
# and the mean-variance relationship is accommodated using precision weights calculated by the voom function.
#
# When operating on a DGEList-object, voom converts raw counts to log-CPM values by automatically extracting 
# library sizes and normalisation factors from x itself. Additional normalisation to log-CPM values can be 
# specified within voom using the normalize.method argument.
#
# The mean-variance relationship of log-CPM values for this dataset is shown in the left-hand panel of the 
# next figure. Typically, the voom-plot shows a decreasing trend between the means and variances resulting 
# from a combination of technical variation in the sequencing experiment and biological variation amongst the 
# replicate samples from different cell populations. Experiments with high biological variation usually result 
# in flatter trends, where variance values plateau at high expression values. Experiments with low biological 
# variation tend to result in sharp decreasing trends.
#
# Moreover, the voom-plot provides a visual check on the level of filtering performed upstream. 
# If filtering of lowly-expressed genes is insufficient, a drop in variance levels can be observed at the low end 
# of the expression scale due to very small counts. If this is observed, one should return to the earlier filtering 
# step and increase the expression threshold applied to the dataset.
#
# Where sample-level variation is evident from earlier inspections of the MDS plot, the  voomWithQualityWeights 
# function can be used to simultaneously incorporate sample-level weights together with the abundance dependent weights 
# estimated by voom (Liu et al. 2015). For an example of this approach, see Liu et al. (2016) (Liu et al. 2016).
#
# NOTE: Be sure to run these one at a time!
#

ar(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

# Figure 4: Means (x-axis) and variances (y-axis) of each gene are plotted to show the dependence between 
# the two befor e voom is applied to the data (left panel) and how the trend is removed after voom precision
# weights are applied to the data (right panel)
#
# The plot on the left is created within the voom function which extracts residual variances from fitting 
# linear models to log-CPM transformed data. Variances are then rescaled to quarter-root variances 
# (or square-root of standard deviations) and plotted against the mean expression of each gene. 
# The means are log2-transformed mean-counts with an offset of 2. The plot on the right is created using 
# plotSA which plots log2 residual standard deviations against mean log-CPM values. The average log2 
# residual standard deviation is marked by a horizontal blue line. 
#
# In both plots, each black dot represents a gene and a red curve is fitted to these points.
#
# Note that the other data frames stored within the DGEList-object that contain gene- and sample-level 
# information, are retained in the EList-object v created by voom. The v$genes data frame is equivalent 
# to x$genes, v$targets is equivalent to x$samples, and the expression values stored in v$E is analogous 
# to x$counts, albeit on a transformed scale. In addition to this, the  voom EList-object has a matrix 
# of precision weights v$weights and stores the design matrix in  v$design.
#

###################################################################################
# 6.3 Fitting linear models for comparisons of interest
###################################################################################

# Linear modelling in limma is carried out using the lmFit and contrasts. Fit functions originally 
# written for application to microarrays. The functions can be used for both microarray and RNA-seq 
# data and fit a separate model to the expression values for each gene. 
#
# Next, empirical Bayes moderation is carried out by borrowing information across all the genes to 
# obtain more precise estimates of gene-wise variability (Smyth 2004). The model’s residual variances 
# are plotted against average expression values in the next figure. It can be seen from this plot that 
# the variance is no longer dependent on the mean expression level.
#

###################################################################################
# 6.4 Examining the number of DE genes
###################################################################################

# For a quick look at differential expression levels, the number of significantly up- and down-regulated 
# genes can be summarised in a table. 
# 
# Significance is defined using an adjusted p-value cutoff that is set at 5% by default. 
#
# For the comparison between expression levels in basal and LP, 4,648 genes are found to be down-regulated 
# in basal relative to LP and 4,863 genes are up-regulated in basal relative to LP – a total of 9,511 DE genes. 
# A total of 9,598 DE genes are found between basal and ML (4,927 down- and 4,671 up-regulated genes), and a 
# total of 5,652 DE genes are found between LP and ML (3,135 down- and 2,517 up-regulated). 
# 
# The larger numbers of DE genes observed for comparisons involving the basal population are consistent with 
# our observations from the MDS plots.
#

summary(decideTests(efit))

# Some studies require more than an adjusted p-value cut-off. For a stricter definition on significance, 
# one may require log-fold-changes (log-FCs) to be above a minimum value. 
# 
# The treat method (McCarthy and Smyth 2009) can be used to calculate p-values from empirical Bayes 
# moderated t-statistics with a minimum log-FC requirement. The number of differentially expressed genes 
# are reduced to a total of 3,648 DE genes for basal versus LP, 3,834 DE genes for basal versus ML, and 
# 414 DE genes for LP versus ML when testing requires genes to have a log-FC that is significantly greater 
# than 1 (equivalent to a 2-fold difference between cell types on the original scale).
#
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

# Genes that are DE in multiple comparisons can be extracted using the results from decideTests, 
# where:
#     0s represent genes that are not DE, 
#     1s represent genes that are up-regulated, and 
#    -1s represent genes that are down-regulated. 
# A total of 2,784 genes are DE in both basal versus LP and basal versus ML, twenty of which are 
# listed below. The write.fit function can be used to extract and write results for all three 
# comparisons to a single output file.
#
de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

head(tfit$genes$SYMBOL[de.common], n=20)

vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

write.fit(tfit, dt, file="results.txt")


###################################################################################
# 6.5 Examining individual DE genes from top to bottom
###################################################################################
#
# The  top DE genes can be listed using topTreat for results using treat (or topTable for results using eBayes). 
# By default topTreat arranges genes from smallest to largest adjusted p-value with associated gene information, 
# log-FC, average log-CPM, moderated t-statistic, raw and adjusted p-value for each gene. The number of top genes 
# displayed can be specified, where  n=Inf includes all genes. Genes Cldn7 and Rasef are amongst the top DE genes 
# for both basal versus LP and basal versus ML.
# 

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)

head(basal.vs.lp)

head(basal.vs.ml)

###################################################################################
# 6.6 Useful graphical representations of differential expression results
###################################################################################

# To summarise results for all genes visually, mean-difference plots, which display 
# og-FCs from the linear model fit against the average log-CPM values can be generated 
# using the plotMD function, with the differentially expressed genes highlighted.
#

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))

# Glimma extends this functionality by providing an interactive mean-difference plot via the
# glMDPlot function. The output of this function is an html page, with summarised results in the 
# left panel (similar to what is output by plotMD), and the log-CPM values from individual samples 
# for a selected gene in the right panel, with a table of results below the plots. This interactive 
# display allows the user to search for particular genes based on the annotation provided (e.g. Gene symbol 
# identifier), which is not possible in a static R plot.

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
        side.main="ENTREZID", counts=lcpm, groups=group, launch=FALSE)

# The mean-difference plot generated by the command above is available online 
# (see http://bioinf.wehi.edu.au/folders/limmaWorkflow/glimma-plots/MD-Plot.html). 
# The interactivity provided by the Glimma package allows additional information to be presented in a 
# single graphical window. Glimma is implemented in R and Javascript, with the R code generating the 
# data which is converted into graphics using the Javascript library D3 (https://d3js.org), with the 
# Bootstrap library handling layouts and Datatables generating the interactive searchable tables. 
# This allows plots to be viewed in any modern browser, which is convenient for including them as linked 
# files from an Rmarkdown report of the analysis.

# Plots shown previously include either all of the genes that are expressed in any one condition 
# (such as the Venn diagram of common DE genes or mean-difference plot) or look at genes individually 
# (log-CPM values shown in right panel of the interactive mean-difference plot). Heatmaps allow users to 
# look at the expression of a subset of genes. This can be give useful insight into the expression of 
# individual groups and samples without losing perspective of the overall study when focusing on individual 
# genes, or losing resolution when examining patterns averaged over thousands of genes at the same time.

# A heatmap is created for the top 100 DE genes (as ranked by adjusted p-value) from the basal versus LP 
# contrast using the heatmap.2 function from the gplots package. The heatmap correctly clusters samples 
# by cell type and reorders the genes into blocks with similar expression patterns. From the heatmap, we observe 
# that the expression of ML and LP samples are very similar for the top 100 DE genes between basal and LP.

library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

# Figure 6: Heatmap of log-CPM values for top 100 genes DE in basal versus LP
# Expression across each gene (or row) have been scaled so that mean expression is zero and standard 
# deviation is one. Samples with relatively high expression of a given gene are marked in red and samples 
# with relatively low expression are marked in blue. 
#
# Lighter shades and white represent genes with intermediate expression levels. 
# Samples and genes have been reordered by the method of hierarchical clustering. 
# A dendrogram is shown for the sample clustering.
#

###################################################################################
# 7 Gene set testing with camera
###################################################################################

# We finish off this analysis with some gene set testing by applying the camera method (Wu and Smyth 2012) 
# to the c2 gene signatures from Broad Institute’s MSigDB c2 collection (Subramanian et al. 2005) 
# that have been adapted for mouse and are available as Rdata objects 
# from http://bioinf.wehi.edu.au/software/MSigDB/. 
# 
# Other useful gene sets derived from MSigDB for both human and mouse, such as the hallmark gene sets, 
# are also available from this site. C2 gene sets have been curated from online databases, publications and 
# domain experts, and hallmark gene sets are selected from MSigDB to have well-defined biological states 
# or processes.

# This part of the script does not nwork



