# Taiji: A Multi-Omic Framework for Gene Regulatory Network Inference

## Overview

**Taiji** is a systems biology pipeline developed to integrate diverse high-throughput omics data—particularly RNA-seq, ATAC-seq, and chromatin interaction data—to infer genome-wide transcription factor (TF) activity and construct gene regulatory networks (GRNs). It provides a quantitative measure of TF influence across regulatory landscapes using personalized PageRank, making it a powerful tool for identifying master regulators in complex biological contexts.

## Key Features

- Integration of **RNA-seq**, **ATAC-seq**, and **3D chromatin structure** data.
- **EpiTensor**-based inference of enhancer–promoter interactions.
- Motif scanning of accessible chromatin regions for TF binding site prediction.
- Personalized **PageRank algorithm** to quantify TF activity across the network.
- Construction of directed gene regulatory networks highlighting TF–target relationships.
- Identifies **differentially active transcription factors** across conditions (e.g., disease vs. control).

## Taiji Pipeline Steps

1. **Chromatin Accessibility Analysis**  
   - ATAC-seq peaks are scanned for TF binding motifs.
   - Peaks are linked to target genes using proximity and EpiTensor-based enhancer–promoter prediction.

2. **Gene Expression Profiling**  
   - RNA-seq quantifies gene expression and provides dynamic transcriptional output.

3. **Network Construction**  
   - A directed regulatory network is built where TFs connect to target genes via motif occurrence in open chromatin.

4. **TF Activity Scoring (PageRank)**  
   - Personalized PageRank is applied to rank TFs based on their global influence in driving gene expression changes.

5. **Differential Analysis**  
   - Compares TF activity between conditions (e.g., RA vs. OA) to highlight disease-specific regulators.

## Application in RA vs. OA FLS

In our study of fibroblast-like synoviocytes (FLS), Taiji revealed:

- **RA-specific TF hubs**: STAT1, STAT3, RELA, and JUNB were highly active in RA FLS, embedded in densely connected TLR pathways.
- **OA as a control**: OA FLS showed activity of TFs linked to tissue maintenance, antioxidant defense, and metabolic homeostasis.
- **Repressed TFs in RA**: Developmental TFs like POU3F3 and NR0B1 were suppressed, suggesting a loss of homeostatic identity.
- **Multi-omic integration**: Taiji linked changes in chromatin accessibility and enhancer–promoter architecture to shifts in transcriptional output and TF activity.

## Resources

- **Taiji Paper**: Zhang K, et al. *Sci Adv.* 2019;5(3):eaav3262. doi:[10.1126/sciadv.aav3262](https://doi.org/10.1126/sciadv.aav3262)
- **EpiTensor**: Computational tool used within Taiji for 3D chromatin interaction modeling.

## Citation

Zhang K, Wang M, Zhao Y, Wang W. *Taiji: System-level identification of key transcription factors reveals transcriptional waves in mouse embryonic development*. Sci Adv. 2019;5(3):eaav3262. doi:10.1126/sciadv.aav3262

