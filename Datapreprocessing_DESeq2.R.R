
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("dplyr")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
query
GDCdownload(query)
brca_data <- GDCprepare(query)
assay(brca_data)[1:5, 1:5]
colData(brca_data)
rowData(brca_data)
table(brca_data$shortLetterCode)
tumor_data <- brca_data[, brca_data$shortLetterCode == "TP"]
normal_data <- brca_data[, brca_data$shortLetterCode == "NT"]
subtypes <- TCGAquery_subtype(tumor = "BRCA")
head(subtypes)
colnames(subtypes)



tumor_data <- brca_data[, brca_data$shortLetterCode == "TP"]
dim(tumor_data)
expr_patients <- substr(colnames(tumor_data), 1, 12)

subtypes <- TCGAquery_subtype(tumor = "BRCA")

subtypes_matched <- subtypes[subtypes$patient %in% expr_patients, ]
nrow(subtypes_matched)

subtypes_ordered <- subtypes_matched[
  match(expr_patients, subtypes_matched$patient),
]
colData(tumor_data)$PAM50 <- subtypes_ordered$BRCA_Subtype_PAM50

table(colData(tumor_data)$PAM50)
saveRDS(tumor_data, file = "TCGA_BRCA_Tumor_STARcounts_PAM50.rds")
saveRDS(subtypes_ordered, file = "TCGA_BRCA_PAM50_subtypes.rds")
tumor_data

getwd()


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
library(DESeq2)
library(ggplot2)

tumor_data <- readRDS("TCGA_BRCA_Tumor_STARcounts_PAM50.rds")
tumor_data
count_matrix <- assay(tumor_data)
col_data <- colData(tumor_data)
col_data$PAM50 <- factor(col_data$PAM50)
keep <- rowSums(count_matrix >= 10) >= 5
count_matrix <- count_matrix[keep, ]
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ PAM50
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeqDataSet(
  tumor_data,
  design = ~ PAM50
)
table(is.na(colData(tumor_data)$PAM50))
tumor_data <- tumor_data[, !is.na(colData(tumor_data)$PAM50)]
table(is.na(colData(tumor_data)$PAM50))
table(colData(tumor_data)$PAM50)
dds <- DESeqDataSet(
  tumor_data,
  design = ~ PAM50
)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
resultsNames(dds)
res_LumA_vs_Basal <- results(
  dds,
  name = "PAM50_LumA_vs_Basal"
)
res_LumA_vs_Basal <- lfcShrink(
  dds,
  coef = "PAM50_LumA_vs_Basal",
  res = res_LumA_vs_Basal,
  type = "apeglm"
)
BiocManager::install("apeglm")
library(apeglm)
res_LumA_vs_Basal <- lfcShrink(
  dds,
  coef = "PAM50_LumA_vs_Basal",
  res = res_LumA_vs_Basal,
  type = "apeglm"
)
summary(res_LumA_vs_Basal)
res_ordered <- res_LumA_vs_Basal[order(res_LumA_vs_Basal$padj), ]
head(res_ordered)
saveRDS(res_LumA_vs_Basal, "DESeq2_LumA_vs_Basal.rds")

write.csv(
  as.data.frame(res_LumA_vs_Basal),
  "DESeq2_LumA_vs_Basal.csv"
)
getwd()
deg_final <- subset(
  as.data.frame(res_LumA_vs_Basal),
  padj < 0.05 & abs(log2FoldChange) > 1
)
write.csv(deg_final, "Final_DEGs_LumA_vs_Basal.csv")
getwd()