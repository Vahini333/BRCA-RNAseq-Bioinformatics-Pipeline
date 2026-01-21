
install.packages("pheatmap")
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
res_LumA_vs_Basal <- readRDS("DESeq2_LumA_vs_Basal.rds")
deg <- res_LumA_vs_Basal[
  which(res_LumA_vs_Basal$padj < 0.05 &
          abs(res_LumA_vs_Basal$log2FoldChange) > 1),
]

tumor_data <- readRDS("TCGA_BRCA_Tumor_STARcounts_PAM50.rds")
list.files()
getwd()
tumor_data
table(colData(tumor_data)$PAM50)
library(DESeq2)

dds <- DESeqDataSet(
  tumor_data,
  design = ~ PAM50
)
table(colData(tumor_data)$PAM50, useNA = "ifany")
tumor_data <- tumor_data[, !is.na(colData(tumor_data)$PAM50)]
table(colData(tumor_data)$PAM50)
colData(tumor_data)$PAM50 <- droplevels(colData(tumor_data)$PAM50)
library(DESeq2)

dds <- DESeqDataSet(
  tumor_data,
  design = ~ PAM50
)
dds <- DESeq(dds)
saveRDS(dds, "dds_TCGA_BRCA.rds")
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "vsd_TCGA_BRCA.rds")
ls()
library(DESeq2)

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, "vsd_TCGA_BRCA.rds")
list.files()
res <- results(
  dds,
  contrast = c("PAM50", "LumA", "Basal")
)
deg <- res[
  which(res$padj < 0.05 &
          abs(res$log2FoldChange) > 1),
]
top_genes <- head(
  rownames(deg[order(deg$padj), ]),
  50
)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))
annotation_col <- data.frame(
  PAM50 = colData(dds)$PAM50
)

rownames(annotation_col) <- colnames(mat_scaled)
library(pheatmap)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Heatmap of Top 50 DEGs: LumA vs Basal"
)
png("TCGA_BRCA_LumA_vs_Basal_Heatmap.png", width = 1000, height = 1200)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Heatmap of Top 50 DEGs: LumA vs Basal"
)
dev.off()
getwd()
png("TCGA_BRCA_LumA_vs_Basal_Heatmap.png",
    width = 1200,
    height = 1400,
    res = 300)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_method = "complete",
  main = "Heatmap of Top 50 DEGs: LumA vs Basal"
)

dev.off()
