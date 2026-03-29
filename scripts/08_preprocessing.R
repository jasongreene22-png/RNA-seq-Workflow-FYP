## Read in library packages required for differntial expression analysis

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


if (!require("apeglm")) {
  BiocManager::install("apeglm")
}
library(apeglm)


if (!require("DESeq2")) {
  BiocManager::install("DESeq2")
}
library(DESeq2)


if (!require("tidyverse")) {
  install.packages("tidyverse")
}
library(tidyverse)


if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)


if (!require("pheatmap")) {
  install.packages("pheatmap")
}
library(pheatmap)


if (!require("ggrepel")) {
  install.packages("ggrepel")
}
library(ggrepel)


if (!require("msigdbr")) {
  install.packages("msigdbr")
}
library(msigdbr)

if (!require("tibble")) {
  install.packages("tibble")
}
library(tibble)

if (!require("fgsea")) {
  BiocManager::install("fgsea")
}
library(fgsea)

if (!require("devtools")) {
  BiocManager::install("devtools")
}
library(devtools)

library(DESeq2)


## Read in Feature Counts text file 
study1counts <- read.table("Study1_featurecounts.txt", header = TRUE)
row.names(study1counts) <- study1counts$Geneid
View(study1counts)



## Remove column names 
study1counts <- study1counts[ , -c(1:6) ]
head(study1counts)

#Repace column names using gsub and assign samples to groups
orig_names <- names(study1counts)
orig_names

names(study1counts) <- gsub("(WT_rep|R47H_rep)([1-9]).*" , "\\1\\2" , orig_names)
names(study1counts)

metadata_study1 <- data.frame(Week = gsub(".*_(week[0-9]).*", "\\1", names(study1counts)), Genotype = gsub(".*_(WT|R47H)_rep[0-9]", "\\1", names(study1counts)),row.names = names(study1counts))
head(metadata_study1)
View(metadata_study1)


## Ensure sample names match in both files
all(colnames(study1counts) %in% rownames(metadata_study1))
identical(colnames(study1counts), rownames(metadata_study1))



