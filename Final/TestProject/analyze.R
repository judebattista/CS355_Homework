source("https://bioconductor.org/biocLite.R")
biocLite("Rsubread")
biocLite("limma")
biocLite("Glimma")
biocLite("edgeR")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
biocLite("KEGGREST")
biocLite("Rsubread")
biocLite("biomaRt")
biocLite("RColorBrewer")
biocLite("gplots")
biocLite("FGNet")
biocLite("CoRegNet")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
biocLite("BSgenome.Celegans.UCSC.ce11")

library("limma")
library("Glimma")
library("edgeR")
library("pathview")
library("gage")
library("gageData")
library("KEGGREST")
library("Rsubread")
library("biomaRt")
library("RColorBrewer")
library("gplots")
library(gage)
library(FGNet)
library("BSgenome.Celegans.UCSC.ce11")
library(AnnotationDbi)
library(CoRegNet)

# next step was to download celegans dna from the the internet then build an index for it
buildindex( basename="celegans_index", reference="./caenorhabditis_elegans.PRJNA13758.WBPS13.genomic.fa" )

#Next I made a Targets.txt file with the list of files to process:
fileConn <- file("Targets.txt")
writeLines(c("CellType\tFilename1\tFilename2\tOutputFile",
             
             "FRK1\t../PutzkeData/frk_1_sample_1/concatenated_fastqs/Sample_mutant-1.concatenated.R1.fastq.gz\t../PutzkeData/frk_1_sample_1/concatenated_fastqs/Sample_mutant-1.concatenated.R2.fastq.gz\tfrk1.bam",
             
             "FRK2\t../PutzkeData/frk_1_sample_2/fastq/frk1-2_ATTCAGAA-GGCTCTGA_AC9FN5ANXX_L008_001.R1.fastq.gz\t../PutzkeData/frk_1_sample_2/fastq/frk1-2_ATTCAGAA-GGCTCTGA_AC9FN5ANXX_L008_001.R2.fastq.gz\tfrk2.bam",
             
             "FRK3\t../PutzkeData/frk_1_sample_3/concatenated_fastqs/Sample_mutant-2.concatenated.R1.fastq.gz\t../PutzkeData/frk_1_sample_3/concatenated_fastqs/Sample_mutant-2.concatenated.R2.fastq.gz\tfrk3.bam",
             
             "WT1\t../PutzkeData/WT_sample_1/concatenated_fastqs/Sample_wildtype-1.concatenated.R1.fastq.gz\t../PutzkeData/WT_sample_1/concatenated_fastqs/Sample_wildtype-1.concatenated.R2.fastq.gz\twt1.bam",
             
             "WT2\t../PutzkeData/WT_sample_2/fastq/wildtype-2_TCTCGCGC-GGCTCTGA_AC9FN5ANXX_L008_001.R1.fastq.gz\t../PutzkeData/WT_sample_2/fastq/wildtype-2_TCTCGCGC-GGCTCTGA_AC9FN5ANXX_L008_001.R2.fastq.gz\twt2.bam",
             
             "WT3\t../PutzkeData/WT_sample_3/concatenated_fastqs/Sample_wildtype-3.concatenated.R1.fastq.gz\t../PutzkeData/WT_sample_3/concatenated_fastqs/Sample_wildtype-3.concatenated.R2.fastq.gz\twt3.bam"
             
), fileConn)

close(fileConn)

# Read in the targets file
targets <- readTargets()

# ran the alignment procedure for all the 6 files
align(index="celegans_index", readfile1=targets$Filename1, readfile2=targets$Filename2, input_format="gzFASTQ", output_format="BAM", output_file=targets$OutputFile, unique=TRUE, indels=5 )

# Downloaded annotation file for c elegans from ftp://ftp.ensemblgenomes.org/pub/metazoa/release-36/gff3/caenorhabditis_elegans

# converted from gff3 to gtf format using the gffread package

# Following command line run from terminal and NOT R to convert gff3 format to gtf format
./gffread Caenorhabditis_elegans.WBcel235.36.gff3 -T -o Caenorhabditis_elegans.WBcel235.36.gtf

# ran feature counts
fc <- featureCounts(files=targets$OutputFile, isPairedEnd=TRUE,annot.ext="Caenorhabditis_elegans.WBcel235.36.gtf", isGTFAnnotationFile=TRUE)

# Read feature counts into a DGEList object using the edgeR package function
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

samplenames <- c("frk1", "frk2", "frk3", "wt1", "wt2", "wt3" )
colnames(x) <- samplenames
group <- as.factor(c("FRK", "FRK", "FRK", "WT", "WT", "WT"))
x$samples$group <- group

lane <- as.factor( c("L002", "L008", "L002", "L002", "L008", "L002") )
x$samples$lane <- lane