#!/bin/bash

# Run FastQC on all gzipped FASTQ files
# Output goes to FastQC directory

mkdir -p FastQC
cd~/FastQC
fastqc *.fastq.gz --outdir FastQC
