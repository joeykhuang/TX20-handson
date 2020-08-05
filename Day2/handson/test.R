if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tximportData")

library("tximport")
library("readr")
library("tximportData")

## read data from tximportData package
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)

## specify the condition, this A,B can be replaced by other term
samples$condition <- factor(rep(c("A","B"),each=3)) 
rownames(samples) <- samples$run 
samples[,c("pop","center","run","condition")]
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz") 
names(files) <- samples$run
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds<-ddsTxi
## filter genes without enough coverage, for example total read depth in all samples are smaller than 10
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

dds <- DESeq(dds) 
res <- results(dds) 
resOrdered <- res[order(res$pvalue),]

summary(res)
write.csv(as.data.frame(resOrdered), file="condition_treated_results.csv")


## Find a biological process you want to explore and get the differentially expressed genes
## and also 

## For example: GSE92572 provide you with count matrix(mapped)
## for example: 