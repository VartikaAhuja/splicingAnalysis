# splicingAnalysis
Codes used to align the samples to human reference, obtain gene and isoform expression profiles, and finally filter cancer-specific isoforms

## RNA-seq Analysis Workflow

## Overview
This repository contains scripts for processing RNA-seq data from raw FASTQ files through alignment, quantification, and differential expression analysis. The workflow supports gene-level, transcript-level, and gene-based transcript-level analyses using **HISAT2**, **StringTie**, **bedtools**, and **DESeq2**.

## Prerequisites

### Software
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
- [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) or [Cutadapt](https://cutadapt.readthedocs.io/)  
- [HISAT2](https://daehwankimlab.github.io/hisat2/)  
- [Samtools](http://www.htslib.org/)  
- [StringTie](https://ccb.jhu.edu/software/stringtie/)  
- [Bedtools](https://bedtools.readthedocs.io/)  
- Python (for `prepDE.py`)  
- R with packages: `DESeq2`, `RColorBrewer`, `pheatmap`

### Reference Files
- Human reference genome (e.g., `hg38.fa`)  
- Gene annotation file (e.g., `gencode.v45.basic.annotation.gtf`)  
- BED file of gene coordinates (`hg38_Genes.bed`)  
- Sample files metadata (`samples_info.txt`)  with columns - sample name, condition and others (like site, batch, stage if available)
- Optional: gene lists (`genelist`, `celf6_genes_literature.txt`)  generated using literature
- Transcript-to-gene mapping (`tid_geneid.txt`) generated using gencode.v45.basic.annotation.gtf

## Workflow

### 1. Quality Control & Trimming
- Run **FastQC** on raw FASTQ files.  
- Trim adapters and low-quality bases using **TrimGalore** or **Cutadapt**.

### 2. Alignment
Code 1_Alignment_and_count_generation.sh can be used (after changing sample names)
- Build HISAT2 index from `hg38.fa`.  
- Align paired-end reads with HISAT2 (`--dta` ensures compatibility with StringTie).  
- Convert SAM → BAM, sort, and index using Samtools.

### 3. Quantification
- Use **bedtools multicov** to generate **gene counts**.  
- Run **StringTie** with annotation GTF to generate transcript assemblies and abundance estimates.  
- Use `prepDE.py` to generate **transcript counts**.

### 4. Differential Expression Analysis
Choose one of the following scripts depending on the analysis type:

- **Code 2 (Gene-level DESeq2)**  
  Identify differentially expressed genes, focusing on splicing factors.

- **Code 3 (Transcript-level DESeq2)**  
  Identify tumor-specific isoforms based on prevalence thresholds (e.g., ≥80% tumor, ≤20% normal).

- **Code 4 (Gene-based Transcript DESeq2)**  
  Identify upregulated splicing factors and filter downstream transcripts regulated by them.

## Inputs & Outputs

### Inputs
- FASTQ files  
- Reference genome (`hg38.fa`)  
- Annotation files (`.gtf`, `.bed`)  
- Sample metadata (`samples_info.txt`)  
- Optional gene lists (`genelist`, `celf6_genes_literature.txt`)  

### Outputs
- Gene counts (`gene_counts.txt`)  
- Transcript counts (`transcript_count_matrix.csv`)  
- DESeq2 results (`deseq2_all_results.txt`, filtered subsets)  
- Plots (heatmaps, PCA, MA plots)  
- Tumor-specific isoform tables (`upregulated_isoforms_all.tsv`, `celf6_Dtranscripts_80T_20N.tsv`)  

## Usage
1. Run **Code 1** to generate gene and transcript counts.  
2. Run **one of Codes 2–4** depending on the analysis type.  
3. Ensure sample metadata (`samples_info.txt`) matches column names in count matrices.  
4. Adjust thresholds (`log2FoldChange`, `pvalue`, prevalence cutoffs) as needed.  

## Notes
- Always verify that sample names in `samples_info.txt` match column names in count matrices.  
- Modify prevalence thresholds (`n1`, `n2`) to suit your experimental design.  
- Gene lists can be customized depending on the biological question.