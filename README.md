# BRCA-RNAseq-Bioinformatics-Pipeline
Uncovering Key Differentially Expressed Genes in Breast Cancer Subtypes Using RNA-Seq

Overview:

This project investigates molecular differences between breast cancer subtypes using RNA sequencing (RNA-Seq) data from The Cancer Genome Atlas (TCGA-BRCA). By comparing Luminal A and Basal-like tumors, we identify differentially expressed genes (DEGs), visualize transcriptomic patterns, and perform functional enrichment analyses to uncover subtype-specific biological pathways. The workflow emphasizes reproducibility, statistical rigor, and biological interpretability.

Data Sources:

Raw Data: TCGA-BRCA RNA-Seq STAR counts

Metadata: PAM50 molecular subtype annotations

Source: The Cancer Genome Atlas (TCGA)

Usage:

Step 1: Data Preprocessing

Load raw count data

Filter low-expression genes

Prepare sample metadata

Step 2: Differential Expression Analysis

Construct DESeqDataSet

Normalize counts and estimate dispersions

Identify DEGs (adjusted p-value < 0.05)

Step 3: Visualization

Volcano plot to highlight significant DEGs

Heatmap of top 50 DEGs using variance-stabilized data

Step 4: Functional Enrichment

Gene Ontology (GO) Biological Process analysis

KEGG pathway enrichment analysis

Results:

Differential Expression:

Thousands of genes showed significant differential expression

Clear upregulation in Luminal A and downregulation in Basal-like tumors

Clustering:

Heatmaps demonstrate strong subtype-specific expression patterns

Samples cluster primarily by PAM50 subtype

Functional Enrichment:

Key pathways identified:

Cell cycle regulation,
Hormone signaling,
Cytokine–cytokine receptor interaction,
Calcium signaling,
Epidermal development.



<img width="2400" height="1800" alt="Volcano_LumA_vs_Basal1" src="https://github.com/user-attachments/assets/8298fca6-7452-44f8-9825-5e869c22a0d1" />

                                    Displays DEGs with log2 fold change vs adjusted p-value

<img width="1200" height="1400" alt="TCGA_BRCA_LumA_vs_Basal_Heatmap" src="https://github.com/user-attachments/assets/0bf3bfba-16de-4c60-8246-b8ce3d2295d7" />

                                    Top 50 DEGs clustered by samples and genes

<img width="2400" height="1800" alt="KEGG_Enrichment_LumA_vs_Basal" src="https://github.com/user-attachments/assets/f9c8c29c-3d49-4f17-a30b-3feced7d162b" />

                                            Enriched signaling pathways

<img width="2400" height="1800" alt="GO_BP_Enrichment_LumA_vs_Basal" src="https://github.com/user-attachments/assets/7131c6ac-63b5-40b6-badf-5bf9e016b170" />

                                              Enriched biological processes
