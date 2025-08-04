# Taiji Output Documentation

## Overview
Taiji generates genome-scale transcriptional regulatory networks by integrating chromatin accessibility (ATAC-seq), transcriptomic profiles (RNA-seq), and chromatin interaction data (EpiTensor). The output files include PageRank-based transcription factor (TF) activity scores, TF-target regulatory relationships, and metadata for downstream network analysis.

---

## Key Output Files

### 1. PageRank Scores
- **Description**: Global influence scores for all transcription factors per sample.
- **Columns**:
  - `TF`: Transcription factor gene symbol
  - `PageRank_Sample1`, `PageRank_Sample2`, ... : PageRank scores across samples
- **Use**:
  - Principal Component Analysis (PCA)
  - Differential TF activity analysis
  - Heatmap generation

---

### 2. Network Edges
- **Description**: Directed TF-target interaction edges with associated weights.
- **Columns**:
  - `Source_TF`
  - `Target_Gene`
  - `Edge_Weight`
- **Use**:
  - Network visualization in iGraph
  - Subnetwork extraction (e.g., top 50 TFs)
  - Hub TF identification
  - File size too large to upload on GitHub

---

### 3. Motif Binding Sites
- **Description**: BED file of predicted TF binding sites in accessible regions.
- **Columns**:
  - `chrom`, `start`, `end`, `TF`, `score`, `strand`
- **Use**:
  - IGV track generation
  - Cross-reference with ATAC-seq peaks
  - File size too large to upload on GitHub

---

### 4. Enhancer and Promoter Links
- **Description**: Predicted enhancer-promoter interactions inferred from EpiTensor.
- **Columns**:
  - `Enhancer_Location`
  - `Promoter_Location`
  - `Target_Gene`
  - `Confidence_Score`
- **Use**:
  - Functional annotation of distal regulatory elements
  - Network edge refinement
  - File size too large to upload on GitHub

---

### 5. MetaData
- **Description**: Metadata file containing sample names, groupings (RA vs OA), and file paths used for each analysis.
- **Use**:
  - Downstream grouping and statistical comparisons
  - Batch reproducibility
  - File size too large to upload on GitHub

---

## Applications of Taiji Output
- TF-centric PCA and clustering of RA vs OA samples
- Identification of top differentially active TFs (e.g., STAT3, RELA, STAT1)
- Pathway enrichment using TF target sets
- Inference of regulatory logic in Toll-like receptor and TGF-β signaling
- Visualization of transcriptional hubs and subnetworks using iGraph

---

## Notes
- Ensure all genomic coordinates are in `hg38` for consistency across datasets.
- When comparing PageRank scores between RA and OA, use non-parametric statistical tests (e.g., Wilcoxon) due to non-normal distribution.
- For reproducibility, store each Taiji run in a separate directory labeled by sample ID and condition.

