##RRHO analysis reveals concordance between individual datasets
## Use RRHO2 package
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
## Create heatmap 
RRHO2_heatmap(RRHO_obj)


## See genes that are either up/up or downn/down
RRHO_obj$genelist_uu
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
