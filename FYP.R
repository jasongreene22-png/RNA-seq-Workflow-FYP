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



study1counts <- read.table("Study1_featurecounts.txt", header = TRUE)
row.names(study1counts) <- study1counts$Geneid
View(study1counts)

study4counts <- read.table("featurecounts_Study4.txt", header = TRUE)
row.names(study4counts) <- study4counts$Geneid
View(study4counts)

study1counts <- study1counts[ , -c(1:6) ]
head(study1counts)

study4counts <- study4counts[ , -c(1:6,11)]
head(study4counts)

orig_names <- names(study1counts)
orig_names

names(study1counts) <- gsub("(WT_rep|R47H_rep)([1-9]).*" , "\\1\\2" , orig_names)
names(study1counts)

orig_names4 <- names(study4counts)
orig_names4

names(study4counts) <- gsub("(WT_rep|R47H_rep)([1-9]).*" , "\\1\\2" , orig_names4)
orig_names4 <- names(study4counts)
names(study4counts) <- sub("^trimmed_", "", orig_names4)
names(study4counts)
head(study4counts)


metadata_study1 <- data.frame(Week = gsub(".*_(week[0-9]).*", "\\1", names(study1counts)), Genotype = gsub(".*_(WT|R47H)_rep[0-9]", "\\1", names(study1counts)),row.names = names(study1counts))
head(metadata_study1)
View(metadata_study1)

metadata_study4 <- data.frame (Clone =gsub (".*_(clone[A-D]).*", "\\1", names(study4counts)), Genotype = gsub(".*_(WT|R47H)_rep[0-9]", "\\1", names(study4counts)),row.names = names(study4counts))
head(metadata_study4)
class(study4counts)
View(metadata_study4)

all(colnames(study1counts) %in% rownames(metadata_study1))
identical(colnames(study1counts), rownames(metadata_study1))

identical(colnames(study4counts), rownames(metadata_study4))

## Study1 

dds<- DESeqDataSetFromMatrix(countData = study1counts, colData = metadata_study1, design = ~ Week + Genotype)
colData(dds) %>% head
assay(dds, "counts" ) %>% head

dds$genotype <- factor(dds$Genotype, levels = c("WT","R47H"))
dim(dds)
dim(dds[rowSums(counts(dds)) > 0, ])
dds <- dds[rowSums(counts(dds)) > 0, ]

rld <- rlog(dds, blind=TRUE)
assay(dds) %>% head
assay(rld) %>% head

PCA1<- plotPCA(rld, intgroup= c("Genotype", "Week")) + geom_text(aes(label=name),hjust="inward", vjust=2, size=2)


# Hierarchical clustering plot
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Genotype, rld$Week, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

## Study 4 
metadata_study4$Clone <- droplevels(factor(metadata_study4$Clone))
metadata_study4$Genotype <- droplevels(factor(metadata_study4$Genotype))

metadata_study4 <- metadata_study4[colnames(study4counts), , drop = FALSE]
collapsed_counts <- t(
  rowsum(
    t(study4counts),
    group = metadata_study4$Clone
  )
)

dim(study4counts)      
dim(collapsed_counts)  
colnames(collapsed_counts)

collapsed_coldata <- unique(metadata_study4[, c("Clone", "Genotype")])
rownames(collapsed_coldata) <- collapsed_coldata$Clone
collapsed_coldata <- collapsed_coldata[colnames(collapsed_counts), ]

dds_collapsed <- DESeqDataSetFromMatrix(
  countData = collapsed_counts,
  colData   = collapsed_coldata,
  design    = ~ Genotype
)


rld_4 <- rlog(dds_collapsed, blind=TRUE)
assay(dds_collapsed) %>% head
assay(rld_4) %>% head

PCA2 <-plotPCA(rld_4, intgroup="Genotype", "clone") + geom_text(aes(label=name),hjust="inward", vjust=2, size=2)
      
sampleDists <- dist(t(assay(rld_4)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld_4$Genotype,rld_4$Clone,sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)  


plotCounts(dds, gene="TREM2", intgroup="Genotype")

# Run the differential expression pipeline on the raw counts
dds <- DESeq(dds)
R47HvsWT <-results(dds, contrast=c("Genotype","R47H","WT"), independentFiltering = TRUE , alpha = 0.05)
summary(R47HvsWT)
head(R47HvsWT)

table(R47HvsWT$padj < 0.05)
R47HvsWT[order(R47HvsWT$padj),] %>% head

plotCounts(dds, gene="TREM2", intgroup="genotype") 

write.table(as.data.frame(R47HvsWT[order(R47HvsWT$padj),] ), file="R47H_vs_WT_dge_study1.txt",sep="\t", quote = FALSE)

R47HvsWT_sorted <- R47HvsWT[order(R47HvsWT$padj),]
R47HvsWT_sorted
R47HvsWT_top20 <- head(R47HvsWT_sorted, n=20)
top20_genes_R47HvsWT <- rownames(R47HvsWT_top20)
top20_genes_R47HvsWT

rlog.dge <- rld[top20_genes_R47HvsWT, ] %>% assay

pheatmap(rlog.dge, scale="row", main = "Differential Gene Expression for Study1 (row-based z-score)")

R47HvsWT_tb <- R47HvsWT %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

R47HvsWT_tb <- R47HvsWT_tb %>% mutate(threshold = padj < 0.05)
R47HvsWT_tb <- R47HvsWT_tb %>% arrange(padj) %>% mutate(genelabels = "")
R47HvsWT_tb$genelabels[1:10] <- R47HvsWT_tb $gene[1:10]



volcano1<- ggplot(R47HvsWT_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("R47H vs WT DE gene expression Study 1") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 


R47HvsWT["TIFAB", ]
R47HvsWT["CLU", ]
R47HvsWT["SERPINA1", ]

###### study 4
dds_collapsed <- DESeq(dds_collapsed)
Study4_R47HvsWT <-results(dds_collapsed , contrast=c("Genotype","R47H","WT"), independentFiltering = TRUE , alpha = 0.05)

summary(Study4_R47HvsWT)
head(Study4_R47HvsWT)
table(Study4_R47HvsWT$padj < 0.05)
Study4_R47HvsWT[order(Study4_R47HvsWT$padj),] %>% head


write.table(as.data.frame(Study4_R47HvsWT[order(Study4_R47HvsWT$padj),] ), file="R47H_vs_WT_dge_study4.txt",sep="\t", quote = FALSE)
Study4_R47HvsWT_sorted <- Study4_R47HvsWT[order(Study4_R47HvsWT$padj),]
Study4_R47HvsWT_top20 <- head(Study4_R47HvsWT_sorted, n=20)
Study4_top20_genes_R47HvsWT <- rownames(Study4_R47HvsWT_top20)

Study4_rlog.dge <- rld_4[Study4_top20_genes_R47HvsWT, ] %>% assay
pheatmap(Study4_rlog.dge, scale="row", main = "Differential Gene Expression for Study4 (row-based z-score)")

Study4_R47HvsWT_tb <- Study4_R47HvsWT %>% data.frame() %>% rownames_to_column(var="gene") %>% as_tibble()

Study4_R47HvsWT_tb <- Study4_R47HvsWT_tb %>% mutate(threshold = padj < 0.05)
Study4_R47HvsWT_tb <- Study4_R47HvsWT_tb %>% arrange(padj) %>% mutate(genelabels = "")
Study4_R47HvsWT_tb$genelabels[1:10] <- Study4_R47HvsWT_tb $gene[1:10]

volcano2 <- ggplot(Study4_R47HvsWT_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("R47H vs WT DE gene expression Study 2") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

Study4_R47HvsWT["TREM2", ]
Study4_R47HvsWT["FEZ1", ]
plotCounts(dds_collapsed, gene="TREM2", intgroup="Genotype")


Study4_R47HvsWT["FEZ1", ]


### Overlap Analysis 
Study1_df <- as.data.frame(R47HvsWT)
Study4_df <- as.data.frame(Study4_R47HvsWT)

padj_cutoff <- 0.05



deg1 <- rownames(R47HvsWT)[
  !is.na(R47HvsWT$padj) &
    R47HvsWT$padj < padj_cutoff &
    abs(R47HvsWT$log2FoldChange) > 0
]

length(deg1)

deg2 <- rownames(Study4_R47HvsWT)[
  !is.na(Study4_R47HvsWT$padj) &
    Study4_R47HvsWT$padj < padj_cutoff &
    abs(Study4_R47HvsWT$log2FoldChange) > 0
]

length(deg2)


Sig_both_Studies <- intersect(deg1,deg2)

lfc1 <- Study1_df[Sig_both_Studies, "log2FoldChange"]
lfc2 <- Study4_df[Sig_both_Studies, "log2FoldChange"]

consistent <- (lfc1 > 0 & lfc2 > 0) | (lfc1 < 0 & lfc2 < 0)
Sig_both_Studies <- Sig_both_Studies[consistent]



p1 <- Study1_df[Sig_both_Studies, "pvalue"]
p2 <- Study4_df[Sig_both_Studies, "pvalue"]

valid <- !is.na(p1) & !is.na(p2)

overlap_genes <- Sig_both_Studies[valid]
p1 <- p1[valid]
p2 <- p2[valid]

fisher_p <- mapply(function(x, y) {
  stat <- -2 * (log(x) + log(y))
  pchisq(stat, df = 4, lower.tail = FALSE)
}, p1, p2)


fisher_results <- data.frame(
  Gene = overlap_genes,
  pvalue_Study1 = p1,
  pvalue_Study4 = p2,
  padj_Study1 = Study1_df[overlap_genes, "padj"],
  padj_Study4 = Study4_df[overlap_genes, "padj"],
  log2FC_Study1 = Study1_df[overlap_genes, "log2FoldChange"],
  log2FC_Study4 = Study4_df[overlap_genes, "log2FoldChange"],
  fisher_p = fisher_p
)



fisher_results <- fisher_results[order(fisher_results$fisher_p), ]
top20_genes <- head(fisher_results, 20)

write.table(fisher_results,
                       file = "fisher_results.txt",
                       sep = "\t",
                       row.names = FALSE,
                       quote = FALSE)



library(ggplot2)

plot_df <- subset(fisher_results, Gene %in% c("IRAK3", "CLU"))

# signed combined significance
plot_df$score <- ifelse(
  plot_df$log2FC_Study1 > 0 & plot_df$log2FC_Study4 > 0,
  -log10(plot_df$fisher_p),
  log10(plot_df$fisher_p)
)

# keep order
plot_df$Gene <- factor(plot_df$Gene, levels = c("IRAK3", "CLU"))

plot_df2 <- subset(fisher_results, Gene %in% c("IRAK3", "CLU"))

plot_long <- rbind(
  data.frame(Gene = plot_df2$Gene, Study = "Study 1", log2FC = plot_df2$log2FC_Study1),
  data.frame(Gene = plot_df2$Gene, Study = "Study 2", log2FC = plot_df2$log2FC_Study4)
)

ggplot(plot_long, aes(x = Gene, y = log2FC, fill = Study)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Gene expression of Top 2 DEGs across both studies",
    x = NULL,
    y = "log2 fold change (R47H vs WT)"
  ) +
  theme_minimal()

## Fisher P overlap dataset
study1 <- Study1_df[, c("log2FoldChange", "pvalue")]
study1$gene <- rownames(study1)

study4 <- Study4_df[, c("log2FoldChange", "pvalue")]
study4$gene <- rownames(study4)


merged_df <- merge(
  study1,
  study4,
  by = "gene",
  suffixes = c("_study1", "_study4")
)


merged_df <- merged_df[
  complete.cases(
    merged_df$log2FoldChange_study1,
    merged_df$log2FoldChange_study4,
    merged_df$pvalue_study1,
    merged_df$pvalue_study4
  ),
]

merged_df$same_direction <- sign(merged_df$log2FoldChange_study1) ==
  sign(merged_df$log2FoldChange_study4)

directional_df <- merged_df[merged_df$same_direction, ]


directional_df <- directional_df[
  directional_df$log2FoldChange_study1 != 0 &
    directional_df$log2FoldChange_study4 != 0,
]


directional_df$fisher_p <- apply(
  directional_df[, c("pvalue_study1", "pvalue_study4")],
  1,
  function(p) {
    chisq_stat <- -2 * sum(log(p))
    pchisq(chisq_stat, df = 2 * length(p), lower.tail = FALSE)
  }
)


directional_df$fisher_padj <- p.adjust(directional_df$fisher_p, method = "BH")


directional_df$direction <- ifelse(
  directional_df$log2FoldChange_study1 > 0 & directional_df$log2FoldChange_study4 > 0,
  "up",
  "down"
)


overlap_sig_genes <- directional_df[directional_df$fisher_padj < 0.05, ]


overlap_sig_genes <- overlap_sig_genes[order(overlap_sig_genes$fisher_padj), ]



overlap_up <- overlap_sig_genes[overlap_sig_genes$direction == "up",]
overlap_down <- overlap_sig_genes[overlap_sig_genes$direction == "down", ]


nrow(overlap_sig_genes)
nrow(overlap_up)
nrow(overlap_down)


write.table(overlap_sig_genes,
              file = "Overlap_sig_genes.txt",
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)


##RRHO 
library(dplyr)
library(tibble)
library(RRHO2)

make_rrho_list <- function(res, gene_colname = "Genes") {
  df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(score = -log10(pvalue) * sign(log2FoldChange)) %>%
    filter(!is.na(gene), !is.na(score), is.finite(score)) %>%
    select(gene, score) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  colnames(df) <- c(gene_colname, "DDE")
  as.data.frame(df, stringsAsFactors = FALSE)
}

list1 <- make_rrho_list(R47HvsWT)
list2 <- make_rrho_list(Study4_R47HvsWT)
resultsNames(dds)

list1 <- data.frame(Genes=as.character(list1$Genes), DDE=as.numeric(list1$DDE))
list2 <- data.frame(Genes=as.character(list2$Genes), DDE=as.numeric(list2$DDE))

common_genes <- intersect(list1$Genes, list2$Genes)
list1 <- subset(list1, Genes %in% common_genes)
list2 <- subset(list2, Genes %in% common_genes)

RRHO_obj <- RRHO2_initialize(
  list1, list2,
  labels = c("Study 1 R47HvsWT", "Study 2 R47HvsWT"),
  log10.ind = TRUE,
  method = "hyper"
)

RRHO2_heatmap(RRHO_obj)

uuRRHO_obj$genelist_uu
RRHO_obj$genelist_dd   
RRHO_obj$genelist_ud   
RRHO_obj$genelist_du  

head(RRHO_obj$genelist_uu, 20)

RRHO2_vennDiagram(RRHO_obj, type="uu")
RRHO2_vennDiagram(RRHO_obj, type="dd")


dd_genes <- RRHO_obj$genelist_dd$gene_list_overlap_dd
uu_genes <- RRHO_obj$genelist_uu$gene_list_overlap_uu

dd_genes
length(uu_genes)
length(dd_genes)

dd_table <- data.frame(
  gene       = dd_genes,
  lfc_study1 = res1[dd_genes, "log2FoldChange"],
  padj_study1= res1[dd_genes, "padj"],
  lfc_study4 = res2[dd_genes, "log2FoldChange"],
  padj_study4= res2[dd_genes, "padj"],
  row.names  = NULL
)

head(dd_table)
str(RRHO_obj$genelist_dd, max.level = 2)

uu_table <- data.frame(
  gene       = uu_genes,
  lfc_study1 = res1[uu_genes, "log2FoldChange"],
  padj_study1= res1[uu_genes, "padj"],
  lfc_study4 = res2[uu_genes, "log2FoldChange"],
  padj_study4= res2[uu_genes, "padj"],
  row.names  = NULL
)

head(uu_table)
nrow(uu_table)
nrow(dd_table)


##fgsea  -- link to tutorial followed : https://biostatsquid.com/fgsea-tutorial-gsea/
library(fgsea)
library(readxl)
hallmarks_df <- msigdbr(species = 'Homo sapiens', category = 'H')
hallmark_paths <- split(hallmarks_df$gene_symbol, hallmarks_df$gs_name)

Reactome_df <- msigdbr(species = 'Homo sapiens', category = 'C2', subcollection = 'CP:REACTOME')
reactome_paths <- split(Reactome_df$gene_symbol, Reactome_df$gs_name)

deg_table <- read_excel("mmc1.xlsx",
                        sheet = "Page 10.DEGs_AD") %>%
  transmute(
    gene   = as.character(row.names),
    group  = as.character(groupID),
    coef   = as.numeric(coef),
    fdr    = as.numeric(fdr)
  ) %>%
  filter(!is.na(gene), !is.na(group), !is.na(coef), !is.na(fdr))


rankings1 <- sign(Study1_df$log2FoldChange)*(-log10(Study1_df$pvalue))
names(rankings1) <- rownames(Study1_df)
rankings1

rankings1 <- sort(rankings1, decreasing = TRUE) 
plot(rankings1)
max(rankings1)
min(rankings1)

rankings2 <- sign(Study4_df$log2FoldChange)*(-log10(Study4_df$pvalue))
names(rankings2) <- rownames(Study4_df)
rankings2 <- sort(rankings2, decreasing = TRUE) 
plot(rankings2)
max(rankings2)
min(rankings2)


max_ranking1 <- max(rankings1[is.finite(rankings1)])
min_ranking1 <- min(rankings1[is.finite(rankings1)])
rankings1 <- replace(rankings1, rankings1 > max_ranking1, max_ranking1 * 10)
rankings1 <- replace(rankings1, rankings1 < min_ranking1, min_ranking1 * 10)
rankings1 <- sort(rankings1, decreasing = TRUE)


max_ranking2 <- max(rankings2[is.finite(rankings2)])
min_ranking2 <- min(rankings2[is.finite(rankings2)])
rankings2 <- replace(rankings2, rankings2 > max_ranking2, max_ranking2 * 10)
rankings2 <- replace(rankings2, rankings2 < min_ranking2, min_ranking2 * 10)
rankings2 <- sort(rankings2, decreasing = TRUE)


ggplot(data.frame(gene_symbol = names(rankings1)[1:50], ranks = rankings1[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data.frame(gene_symbol = names(rankings2)[1:50], ranks = rankings2[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


fg1_hall <- fgsea(pathways = hallmark_paths,
                  stats   = rankings1,
                  minSize = 15,
                  scoreType = 'std',
                  maxSize = 500,
                  nproc = 1)

fg2_hall <- fgsea(pathways = hallmark_paths,
                  stats   = rankings2,
                  scoreType = 'std',
                  minSize = 15,
                  maxSize = 500,
                  nproc = 1)


number_of_top_pathways_up <- 10 
topPathwaysUp <- fg1_hall[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg1_hall[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_paths[topPathways], stats = rankings1, fgseaRes = fg1_hall, gseaParam = 0.5)

topPathwaysUp <- fg2_hall[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg2_hall[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(hallmark_paths[topPathways], stats = rankings2, fgseaRes = fg2_hall, gseaParam = 0.5)


fg1_react <- fgsea(pathways = reactome_paths,
                   stats   = rankings1,
                   scoreType = 'std',
                   minSize = 15,
                   maxSize = 500,
                   nproc = 1)

fg2_react <- fgsea(pathways = reactome_paths,
                   stats   = rankings2,
                   scoreType = 'std',
                   minSize = 15,
                   maxSize = 500,
                   nproc = 1)



number_of_top_pathways_up <- 10                       
topPathwaysUp <- fg1_react[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg1_react[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(reactome_paths[topPathways], stats = rankings1, fgseaRes = fg1_react, gseaParam = 0.5)

                    
topPathwaysUp <- fg2_react[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg2_react[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(reactome_paths[topPathways], stats = rankings2, fgseaRes = fg2_react, gseaParam = 0.5)




custom_paths <- deg_table %>%
  filter(fdr < 0.05) %>%
  mutate(set = paste0(group, ifelse(coef > 0, "_AD_up", "_AD_down"))) %>%
  group_by(set) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

custom_paths <- setNames(custom_paths$genes, custom_paths$set)


fg1_custom <- fgsea(pathways = custom_paths, stats = rankings1,
                    scoreType = 'std',
                    minSize = 15,
                    maxSize = 500,
                    nproc = 1)


fg2_custom <- fgsea(pathways = custom_paths, stats = rankings2,
                    scoreType = 'std',
                    minSize = 15,
                    maxSize = 500,
                    nproc = 1)

number_of_top_pathways_up <- 10                       
topPathwaysUp <- fg1_custom[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg1_custom[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(custom_paths[topPathways], stats = rankings1, fgseaRes = fg1_custom, gseaParam = 0.5)

number_of_top_pathways_up <- 10                       
topPathwaysUp <- fg2_custom[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
number_of_top_pathways_down <- 10
topPathwaysDown <- fg2_custom[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(custom_paths[topPathways], stats = rankings2, fgseaRes = fg2_custom, gseaParam = 0.5)




fg1_custom[order(fg1_custom$padj), ][1:10, ]
fg1_customdf <- as.data.frame(fg1_custom)
fg1_custom[fg1_custom$pathway=="MG2.late_AD_down", "leadingEdge"][[1]]


fg2_custom[order(fg2_custom$padj), ][1:10, ]
fg2_customdf <- as.data.frame(fg2_custom)
fg2_custom[fg2_custom$pathway=="MG2.late_AD_up", "leadingEdge"][[1]]


head(fg1_hall[order(fg1_hall$padj), ], 10)
head(fg2_hall[order(fg2_hall$padj), ], 10)

head(fg1_react[order(fg1_react$padj), ], 10)
head(fg2_react[order(fg2_react$padj), ], 10)

sum(fg1_hall[, padj < 0.01])

sum(fg2_hall[, padj < 0.01])


topTerm1 <- fg1_hall[order(padj)][1]$pathway


plotEnrichment(hallmark_paths[[topTerm1]], stats = rankings1) +
  labs(title = paste("Study1:", topTerm1))

plotEnrichment(hallmark_paths[['HALLMARK_INTERFERON_GAMMA_RESPONSE']],
               rankings1) + 
  labs(title = 'HALLMARK pathway: INTERFERON_GAMMA_RESPONSE ') 
 
fg1_hall[fg1_hall$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE", ]



topTerm2 <- fg2_hall[order(padj)][1]$pathway

plotEnrichment(hallmark_paths[[topTerm2]], stats = rankings2) +
  labs(title = paste("Study 2:", topTerm2))

p <- plotEnrichment(
  hallmark_paths[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  rankings2
) + 
  labs(title = "HALLMARK pathway: INTERFERON_GAMMA_RESPONSE") + 
  theme_classic() +
  scale_x_continuous("Rank", breaks = seq(0, 32000, 5000)) +
  scale_y_continuous("Enrichment score (ES)")

p$layers[[1]]$aes_params$colour <- "purple"
p$layers[[1]]$aes_params$linewidth <- 2

p



topTerm2_reactome <- fg2_react [order(padj)][1]$pathway
plotEnrichment(reactome_paths[[topTerm2_reactome]], stats = rankings2) +
  labs(title = paste("Study 2:", topTerm2_reactome))



plotEnrichment(reactome_paths[['REACTOME_INTERFERON_ALPHA_BETA_SIGNALING']],
               rankings1) + 
  labs(title = 'Reactome pathway: INTERFERON_Alpha_B_Signaling') 

fg1_react[fg1_react$pathway == "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"]

fg2_react[fg2_react$pathway == "REACTOME_ECM_PROTEOGLYCANS"]


###Cluster Profiling

library(clusterProfiler)

if (!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
if (!require("ReactomePA")) {
  BiocManager::install("ReactomePA")
}
library(ReactomePA)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db) 
library(dplyr)
library(ggplot2)

ad_genes<- readLines("AD_loci.txt")
ego_ad <- enrichGO(
  gene      = ad_genes,
  OrgDb     = org.Hs.eg.db,
  keyType   = "SYMBOL",     
  ont       = "BP",
  readable  = TRUE
)
s_ego_ad<-clusterProfiler::simplify(ego_ad)
dotplot(s_ego_ad,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of AD related genes")

s_ego_ad %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of AD related genes")

cnetplot(ego_ad,
         showCategory = 5,
         node_label = "category",
         color_category="purple")


#### Overlapped Genes Overrepresentation 
universe <- intersect(rownames(Study1_df), rownames(Study4_df))

ego_overlap <- enrichGO(
  gene      = overlap_sig_genes$gene,
  OrgDb     = org.Hs.eg.db,
  keyType   = "SYMBOL",     
  ont       = "BP",
  universe = universe,
  readable  = TRUE
)

s_ego_overlap<-clusterProfiler::simplify(ego_overlap)
dotplot(ego_overlap,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of overlap DE genes")

s_ego_overlap %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of overlap DE genes")



cnetplot(ego_overlap,
         node_label = "category",
         showCategory = 5 ,
         color_category = "purple")


universe_df<- as.data.frame(universe)

ids<-bitr(overlap_sig_genes$gene, fromType="SYMBOL", 
          toType="ENTREZID", 
          OrgDb="org.Hs.eg.db")



kegg_universe <- bitr(universe, fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")
universe_entrez <- as.character(kegg_universe$ENTREZID)

kgo <- enrichKEGG(gene = ids$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff  = 0.05,
                  universe = universe_entrez
)

mlist<-list(s_ego_overlap,kgo)
names(mlist)<-c("GO-enrich","KEGG-enrich")
mresult<-merge_result(mlist)
dotplot(mresult,showCategory=10)


ego_overlap_up <- enrichGO(gene= overlap_up$gene,
                           OrgDb= org.Hs.eg.db,
                           keyType = "SYMBOL",
                           universe = universe,
                           ont= "BP")

s_ego_overlap_up<-clusterProfiler::simplify(ego_overlap_up)
dotplot(ego_overlap_up,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of overlap up DE genes")

s_ego_overlap_up %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of overlap up DE genes")

nrow(ego_overlap_up)
cnetplot(ego_overlap_up,
         node_label = "all",
         showCategory = 1,
         color_category = "purple")



ego_overlap_down <- enrichGO(gene= overlap_down$gene,
                             OrgDb= org.Hs.eg.db,
                             keyType = "SYMBOL",
                             universe = universe,
                             ont= "BP")


s_ego_overlap_down<-clusterProfiler::simplify(ego_overlap_down)
dotplot(ego_overlap_down,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of overlap down DE genes")


s_ego_overlap_down %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of overlap down DE genes")

cnetplot(ego_overlap_down,
         node_label = "all",
         color_category = "purple")




ad_genes <- unique(ad_genes)
universe <- unique(universe)
overlap_sig_genes <- overlap_sig_genes[!duplicated(overlap_sig_genes$gene), , drop = FALSE]
overlap_ad_dds<- intersect(ad_genes , overlap_sig_genes$gene)
k <- length(overlap_ad_dds)         
n <- nrow(overlap_sig_genes)             
K <- length(ad_genes)             
N <- length(universe )                  
phyper(k-1, K, N-K, n, lower.tail = FALSE)

expected <- (K/N) * n
expected
fold_enrichment <- k / expected
fold_enrichment


AD_overlap_genes<-overlap_sig_genes[overlap_sig_genes$gene %in% ad_genes, 
                                                         c("gene", "fisher_padj","direction")]
write.table(AD_overlap_genes,
            file = "AD_overlap_genes.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)












install.packages("gridExtra")
library(gridExtra)
gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2)
gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2, widths = c(1,1))

grid.arrange(volcano1, volcano2, ncol = 2)
grid.arrange(PCA1, PCA2, ncol = 2)
grid.arrange(fgsea1, fgsea2, ncol = 2)

gridExtra::grid.arrange(up$gtable, down$gtable, ncol = 2)
gridExtra::grid.arrange(h1$gtable, h2$gtable, ncol = 2)





pw <- "MG0.late_AD_down"
le <- fg1_custom$leadingEdge[fg1_custom$pathway == pw][[1]]

le_in <- intersect(le, names(rankings1))
mean(rankings1[le_in] > 0)
length(le_in)

plotEnrichment(custom_paths[[pw]], rankings1)

R47HvsWT["TNF", ]
Study4_R47HvsWT["TNF", ]



##heatmap
topN<-20
r1 <- Study1_df
r1$gene_id <- rownames(r1)

r2 <- Study4_df
r2$gene_id <- rownames(r2)


if ("gene" %in% colnames(overlap_sig_genes)) {
  overlap_sig_genes$gene_id <- overlap_sig_genes$gene
}


Top20DEGsFisher <- head(overlap_sig_genes, topN)


write.table(
  Top20DEGsFisher,
  file = "Top20DEGsFisher.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


top_genes <- Top20DEGsFisher$gene_id


top_genes <- intersect(top_genes, rownames(assay(rld)))
top_genes <- intersect(top_genes, rownames(assay(rld_4)))


mat1 <- assay(rld)[top_genes, , drop = FALSE]
mat2 <- assay(rld_4)[top_genes, , drop = FALSE]


mat1_scaled <- t(scale(t(mat1)))
colnames(mat1_scaled) <- c("WT_W7","WT_W8","WT_W9","R47H_W7","R47H_W8","R47H_W9")
mat2_scaled <- t(scale(t(mat2)))


mat1_scaled[is.na(mat1_scaled)] <- 0
mat2_scaled[is.na(mat2_scaled)] <- 0


tmp <- pheatmap(mat1_scaled, silent = TRUE)
gene_order <- rownames(mat1_scaled)[tmp$tree_row$order]

mat1_scaled <- mat1_scaled[gene_order, , drop = FALSE]
mat2_scaled <- mat2_scaled[gene_order, , drop = FALSE]


r47h_1 <- grep("R47H", colnames(mat1_scaled), value = TRUE)
wt_1   <- setdiff(colnames(mat1_scaled), r47h_1)

mat1_scaled <- mat1_scaled[, c(wt_1, r47h_1), drop = FALSE]
gap1 <- length(wt_1)


geno_map2 <- c(
  cloneA = "WT",
  cloneB = "WT",
  cloneC = "R47H",
  cloneD = "R47H"
)

geno_map2 <- geno_map2[colnames(mat2_scaled)]

ord2 <- c(
  names(geno_map2[geno_map2 == "WT"]),
  names(geno_map2[geno_map2 == "R47H"])
)

mat2_scaled <- mat2_scaled[, ord2, drop = FALSE]
gap2 <- sum(geno_map2[ord2] == "WT")


p1 <- pheatmap(
  mat1_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  gaps_col = gap1,
  main = "Top 20 shared DEGs (same direction; Fisher combined) - Study 1"
)

p2 <- pheatmap(
  mat2_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  gaps_col = gap2,
  main = "Top 20 shared DEGs (same direction; Fisher combined) - Study 2"
)




## KO Overlap#
overlap_sig_genes<- read.table(file = "overlap_sig_genes.txt", header = FALSE, sep = "", quote = "\"'")
sig_genes_sorted<- read.table(file = "sig_genes_sorted.txt", header = FALSE, sep = "", quote = "\"'")
KO_sig_genes<- as.data.frame(sig_genes_sorted)
library(dplyr)
colnames(KO_sig_genes) <- as.character(KO_sig_genes[1,])
colnames(overlap_sig_genes) <- as.character(overlap_sig_genes[1,])

KO_R47H_overlap <- inner_join(KO_sig_genes, overlap_sig_genes, by = "gene")


as.data.frame(KO_R47H_overlap)
nrow(KO_R47H_overlap)

write.table(KO_R47H_overlap,
            file = "KO_R47H_overlap.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

ego_KO_R47H <- enrichGO(gene= KO_R47H_overlap$gene,
                             OrgDb= org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont= "BP")

s_ego_KO_R47H<-clusterProfiler::simplify(ego_KO_R47H)
dotplot(s_ego_KO_R47H,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of overlap of KO and R47H DE genes")

s_ego_KO_R47H %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of overlap of KO and R47H DE genes DE genes")
  





## KO Overlap#
overlap_sig_genes<- read.table(file = "overlap_sig_genes.txt", header = TRUE, sep = "", quote = "\"'")
sig_genes_sorted<- read.table(file = "sig_genes_sorted.txt", header = TRUE, sep = "", quote = "\"'")
KO_sig_genes<- as.data.frame(sig_genes_sorted)
library(dplyr)


KO_R47H_overlap <- inner_join(KO_sig_genes, overlap_sig_genes, by = "gene")


as.data.frame(KO_R47H_overlap)
nrow(KO_R47H_overlap)

write.table(KO_R47H_overlap,
            file = "KO_R47H_overlap.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

ego_KO_R47H <- enrichGO(gene= KO_R47H_overlap$gene,
                             OrgDb= org.Hs.eg.db,
                             keyType = "SYMBOL",
                             ont= "BP")

s_ego_KO_R47H<-clusterProfiler::simplify(ego_KO_R47H)
dotplot(s_ego_KO_R47H,showCategory=20,font.size=10,label_format=70)+
  scale_size_continuous(range=c(1, 7))+
  theme_minimal() +
  ggtitle("GO Enrichment of overlap of KO and R47H DE genes")

s_ego_KO_R47H %>% filter(p.adjust < 0.03) %>%
  ggplot(showCategory = 20,
         aes(GeneRatio, forcats::fct_reorder(Description, GeneRatio))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(1, 7)) +
  theme_minimal() + 
  xlab("Gene Ratio") +
  ylab(NULL) + 
  ggtitle("GO Enrichment of overlap of KO and R47H DE genes DE genes")


library(UpSetR)
library(dplyr)

library(UpSetR)
library(dplyr)


upset_data <- KO_R47H_overlap %>%
  mutate(
    Study_1 = as.integer(!is.na(pvalue_1) & pvalue_1 < 0.05),
    Study_2 = as.integer(!is.na(pvalue_2) & pvalue_2 < 0.05),
    Study_3 = as.integer(!is.na(pvalue_3) & pvalue_3 < 0.05),
    Study_4 = as.integer(!is.na(pvalue_study1) & pvalue_study1 < 0.05),
    Study_5 = as.integer(!is.na(pvalue_study4) & pvalue_study4 < 0.05)
  ) %>%
  dplyr::select(Study_1, Study_2, Study_3, Study_4, Study_5)


upset(as.data.frame(upset_data), 
      sets = c("Study_1", "Study_2", "Study_3", "Study_4", "Study_5"),
      order.by = "freq",
      main.bar.color = "steelblue")


colSums(upset_data)









library(pheatmap)

# -----------------------------
# Study 1: extract and scale
# -----------------------------
mat_1 <- assay(rld)[top20_gene_list, , drop = FALSE]
mat_1_scaled <- t(scale(t(mat_1)))

sample_info_1 <- data.frame(
  sample_original = colnames(mat_1_scaled),
  genotype = ifelse(grepl("R47H", colnames(mat_1_scaled)), "R47H", "WT"),
  week_num = as.numeric(sub(".*week([0-9]+).*", "\\1", colnames(mat_1_scaled))),
  stringsAsFactors = FALSE
)

sample_info_1$sample_label <- paste0(sample_info_1$genotype, "_W", sample_info_1$week_num)
colnames(mat_1_scaled) <- sample_info_1$sample_label

tmp <- pheatmap(mat_1_scaled, silent = TRUE)
gene_order <- rownames(mat_1_scaled)[tmp$tree_row$order]

mat_1_scaled <- mat_1_scaled[gene_order, , drop = FALSE]

sample_info_1 <- sample_info_1[order(sample_info_1$genotype != "WT", sample_info_1$week_num), ]
col_order_1 <- sample_info_1$sample_label
mat_1_scaled <- mat_1_scaled[, col_order_1, drop = FALSE]

gap_1 <- sum(sample_info_1$genotype == "WT")


mat_2 <- assay(rld_4)[top20_gene_list, , drop = FALSE]
mat_2_scaled <- t(scale(t(mat_2)))
mat_2_scaled <- mat_2_scaled[gene_order, , drop = FALSE]

geno_map_2 <- c(
  cloneA = "WT",
  cloneB = "WT",
  cloneC = "R47H",
  cloneD = "R47H"
)

sample_info_2 <- data.frame(
  sample_original = colnames(mat_2_scaled),
  genotype = unname(geno_map_2[colnames(mat_2_scaled)]),
  stringsAsFactors = FALSE
)

sample_info_2 <- sample_info_2[order(sample_info_2$genotype != "WT", sample_info_2$sample_original), ]
sample_info_2$sample_label <- c("WT_1", "WT_2", "R47H_1", "R47H_2")

colnames(mat_2_scaled) <- sample_info_2$sample_label
mat_2_scaled <- mat_2_scaled[, sample_info_2$sample_label, drop = FALSE]

gap_2 <- sum(sample_info_2$genotype == "WT")



h1 <- pheatmap(
  mat_1_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gap_1,
  main = "Top 20 DEGs Significant across both studies - Study 1"
)

h2 <- pheatmap(
  mat_2_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  gaps_col = gap_2,
  main = "Top 20 DEGs Significant across both studies - Study 2"
)









