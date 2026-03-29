## Run Differnetial Expression Analysis on Individual Study
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

## Plot PCA 
plotPCA(rld, intgroup= c("Genotype", "Week")) + geom_text(aes(label=name),hjust="inward", vjust=2, size=2)


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

# Run the differential expression pipeline on the raw counts
dds <- DESeq(dds)
R47HvsWT <-results(dds, contrast=c("Genotype","R47H","WT"), independentFiltering = TRUE , alpha = 0.05)
summary(R47HvsWT)
head(R47HvsWT)

table(R47HvsWT$padj < 0.05)
R47HvsWT[order(R47HvsWT$padj),] %>% head

plotCounts(dds, gene="TREM2", intgroup="genotype") 

write.table(as.data.frame(R47HvsWT[order(R47HvsWT$padj),] ), file="R47H_vs_WT_dge_study1.txt",sep="\t", quote = FALSE)


# Plot Heatmap of the Top 20 DE genes 
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

