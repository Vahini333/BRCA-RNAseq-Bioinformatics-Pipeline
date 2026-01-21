res <- readRDS("DESeq2_LumA_vs_Basal.rds")
res_df <- as.data.frame(res)
deg <- subset(
  res_df,
  padj < 0.05 & abs(log2FoldChange) > 1
)
nrow(deg)
gene_symbols <- rownames(deg)
head(gene_symbols)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
gene_entrez <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(gene_symbols)
gene_symbols_clean <- gsub("\\..*", "", gene_symbols)
head(gene_symbols_clean)
gene_entrez <- bitr(
  gene_symbols_clean,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
gene_entrez_unique <- unique(gene_entrez$ENTREZID)
length(gene_entrez_unique)
kegg_results <- enrichKEGG(
  gene         = gene_entrez_unique,
  organism     = "hsa",
  pvalueCutoff = 0.05
)
head(as.data.frame(kegg_results))
dotplot(kegg_results, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment: Luminal A vs Basal-like")
ggsave("KEGG_Enrichment_LumA_vs_Basal.png", width = 8, height = 6, dpi = 300)



go_results <- enrichGO(
  gene          = gene_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)
dotplot(go_results, showCategory = 10) +
  ggtitle("GO Biological Processes Enriched in DEGs")
ggsave("GO_BP_Enrichment_LumA_vs_Basal.png", width = 8, height = 6, dpi = 300)
