# RNA-seq and Differential Gene Expression Analysis Basic Workflow

## Purpose
This repository contains template scripts outlining a typical RNA-seq
differential expression analysis workflow.

The scripts are intended as **examples and starting points** for users to
adapt to their own data, computing environment, and analysis goals.


## Workflow overview
The workflow follows a standard RNA-seq analysis structure:

1. Quality control of raw sequencing reads (FastQC,MultiQC,Cutadapt)
2. Read alignment or quantification (STAR,samtools)
3. Generation of gene-level count matrices (GeneBodyCoverage, featureCounts, readistribution)
4. Differential expression analysis (DESeq2)
5. Rank-Rank Hypergeometirc Overlap (RRHO)
6. Fast Gene Set Enrichment Analysis (fgsea)
7. Overrepresentation Analysis (ORA)

## Repository structure
Scripts are organized by analysis step:

- `scripts/`
  - `01__02_qc.sh`  
    Example quality control script (e.g. FastQC, MultiQC)

  - `03_cutadapt_sh`  
    Example Trimming for illumina universal adapters and polyA

  - `04_STAR.sh`  
    Example script for alignment using STAR 

  - `05_Samtools.sh`  
    How to interpret alignemnt data using Samtools

   -  `06_RSeqc.sh`
     Examples of code for read distrubtion, geneBodyCoverage and inferring strandness
   -  `07_featurecounts.sh`
      Example feature count scripts for forward and reverse stranded libraries

## Useful links for executing commands and understanding output 
- FastQC - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/INSTALL.txt
- MultiQC - https://multiqc.info/
- Cutadapt - https://cutadapt.readthedocs.io/en/stable/guide.html
- STAR - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- Samtools - http://www.htslib.org/doc/samtools.html
- RSeQC - http://rseqc.sourceforge.net/
- featurecounts - http://subread.sourceforge.net/featureCounts.html


## How to use this repository
Users are expected to:
- Read through each script to understand the workflow logic
- Modify file paths, parameters, and tool versions
- Adapt the scripts to their own RNA-seq data and computing environment
- Integrate the scripts into their own pipeline or workflow manager if desired

## Notes
- No test data is included
- File paths and software calls are system-specific
- Scripts are provided for instructional and reference purposes




