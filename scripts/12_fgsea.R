##fgsea  -- link to tutorial followed : https://biostatsquid.com/fgsea-tutorial-gsea/
library(fgsea)

##Download Hallmark and Reactome Pathways from Msigdb (and others of interest if wanted0
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

## Rank gene lists
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

