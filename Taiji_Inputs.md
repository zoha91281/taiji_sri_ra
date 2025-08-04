# Taiji Input Requirements

## 1. RNA-seq Gene Expression  
- Tab-separated file (TSV)  
- Rows: gene symbols  
- Columns: expression values per sample (e.g., TPM or log2(TPM+1))  
- Header with gene names and sample IDs  

## 2. ATAC-seq Peaks  
- BED file with open chromatin regions  
- Columns: chromosome, start, end  
- Peaks should be merged across all samples to form a consensus peak set  

## 3. Reference Genome Annotation  
- RefGene file or GTF converted to Taiji-compatible format  
- Must include gene name, chromosome, strand, and transcription start site  

## 4. Transcription Factor Motif Database  
- In MEME or JASPAR format  
- Used to scan ATAC-seq peaks for TF binding sites  

## 5. Chromatin Interaction Data (Optional)  
- File linking enhancers to promoters  
- Can be predicted using tools like EpiTensor or derived from Hi-C data  

## 6. Sample Metadata (Optional)  
- Tab-separated file describing each sample’s condition and type  
- Recommended for grouping samples during analysis  

## Preprocessing Requirements  
- Gene expression must be normalized  
- Peaks must be coordinate-sorted  
- Matched RNA-seq and ATAC-seq samples preferred for regulatory modeling  
- All input files must be formatted consistently and error-free  
