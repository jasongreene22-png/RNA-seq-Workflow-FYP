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
