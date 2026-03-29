## Run Overrepresentation analysis on DEG list
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
