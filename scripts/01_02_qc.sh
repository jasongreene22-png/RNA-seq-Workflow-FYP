#!/bin/bash

# Run FastQC on all gzipped FASTQ files
#Use threads (t -) ,where t =  number of CPU threads (i.e. how fast you want to run!) 
#Output goes to FastQC directory

mkdir -p FastQC
cd ~/FastQC
fastqc -t 4 *.fastq.gz --outdir FastQC


#To assess and interpret your FastQC results use this link to help understand  
#https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/07_qc_fastqc_assessment.html


#MultiQC gives will provide easy to interpret data about your raw sequences by giving reports on collection of sample files
# to run multiqc on a directory containing fastqc reports for example run 
mutliqc ~/FastQc

#For more info on how to interpret MultiQC reports visit https://multiqc.info/ 
#To find out what version of multiqc you are using 
multiqc --version 
