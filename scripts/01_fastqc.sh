#!/bin/bash

# Run FastQC on all gzipped FASTQ files
# Output goes to FastQC directory

mkdir -p FastQC
cd~/FastQC
fastqc *.fastq.gz --outdir FastQC

#To run FastQC on multiple files in a directroy use threads (t -) ,where t = number of FastQ  files you want to run
fastqc -t 4 *.fastq.gz --outdir 

#To assess and interpret your FastQC results use this link to help understand  
#https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/07_qc_fastqc_assessment.html


