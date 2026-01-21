library(ggplot2)
res <- readRDS("DESeq2_LumA_vs_Basal.rds")
res_df$Significance <- "Not Significant"

res_df <- as.data.frame(res)

# Remove rows with NA adjusted p-values
res_df <- res_df[!is.na(res_df$padj), ]

# Add gene names as a column
res_df$gene <- rownames(res_df)
res_df$Significance <- "Not Significant"

res_df$Significance[
  res_df$padj < 0.05 & res_df$log2FoldChange > 1
] <- "Upregulated (LumA)"

res_df$Significance[
  res_df$padj < 0.05 & res_df$log2FoldChange < -1
] <- "Downregulated (Basal)"
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated (LumA)" = "red",
    "Downregulated (Basal)" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Luminal A vs Basal-like Breast Cancer",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Gene Status"
  )

volcano_plot
ggsave(
  "Volcano_LumA_vs_Basal.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)
getwd()
