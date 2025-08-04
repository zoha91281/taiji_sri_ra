# Requirements and Tools for Multi-Omic Analysis Using Taiji

## Environments
- **Python**: 3.9
- **R**: 4.2.2
- **Google Colab**: Used for scalable, reproducible workflows

---

## Data Sources
- **Dataset**: NCBI GEO GSE112658
  - RNA-seq (Transcriptomics)
  - ATAC-seq (Chromatin Accessibility)
  - EpiTensor-predicted 3D Chromatin Interactions
  - Reference Genome: hg38

---

## Core Tools and Packages

### 1. **Taiji Pipeline**
- **Version**: 1.2.3
- **Purpose**: Integrates ATAC-seq, RNA-seq, and chromatin loop data
- **Output**: PageRank-based transcription factor influence scores

### 2. **EpiTensor**
- **Purpose**: Predicts enhancer-promoter chromatin interactions

### 3. **scikit-learn**
- **Module**: `sklearn.decomposition.PCA`
- **Version**: 1.6.1
- **Purpose**: Principal Component Analysis on RNA-seq and PageRank matrices

### 4. **GSEAPY**
- **Version**: 1.1.9
- **Purpose**: Gene Set Enrichment Analysis (GSEA) for pathway enrichment

### 5. **python-igraph**
- **Version**: 0.11.0
- **Purpose**: TF-target regulatory network construction and topology analysis

### 6. **cairocffi**
- **Version**: 0.9.1
- **Purpose**: Rendering high-resolution iGraph visualizations

### 7. **pheatmap (R)**
- **Version**: 1.0.13
- **Purpose**: Heatmap visualization of TF PageRank variance

### 8. **ClusterProfiler (R)**
- **Version**: 4.6.0
- **Purpose**: Functional annotation and pathway enrichment visualization

---

## Statistical Methods
- **Wilcoxon Rank-Sum Test**: For identifying differentially active TFs
- **Variance Stabilization & Log Normalization**: Preprocessing RNA-seq
- **Benjamini-Hochberg Correction**: Controlling false discovery rate in GSEA
