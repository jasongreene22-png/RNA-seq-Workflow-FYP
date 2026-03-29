## If performing DE on two indepemdent datasets and wish to see overlap between DEGs 
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

## If interested in genes that are signifciant in EACH individual study
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

## Apply Fisher Combined Probability Test
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



## IF want genes that show overall signficant Fisher P BETWEEN studies
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
