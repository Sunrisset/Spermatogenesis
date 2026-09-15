#!/bin/bash
# ============================================================
# Script name: FASTQs_analysis.sh
# Subject: Analysis of public scRNA-Seq raw data of cells at different stages during mouse spermatogenesis/of germ cells during mouse spermatogenesis
# Purpose: Manage FASTQ files
# Author: Xiaoya
# Date: 2026-09-10
# ============================================================


# ============== Download data. ==============

pwd
# /home/data/SpermatogenicCells/

# ---Download scRNA-seq data.---
# Download scRNA-seq data. It is raw data here. Accession code E-MTAB-6946.  
# Download all fasq files of ten time serials samples.

curl https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6946/ -o sperm_html
# cat sperm_html
grep ".gz" sperm_html > sperm.txt

cat sperm.txt | awk -v FS='"' '{print $8}' > sperm2.txt

awk '{print "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6946/" $0}' sperm2.txt > all_fasq_files.txt

wget -b -i all_fasq_files.txt


# ---Download reference data.---
# In the paper, GRCm38 for mouse + the sequence of the human chromosome 21 (taken from GRCh38).
# For P5-P35, the genetic background of the samples is Mus musculus only. So will use GRCm38.
# For Tc1 and Tc0, the genetic background of the samples is Mus musculus + human chromosome 21. So will use GRCm38 and chro 12 of GRCh38.

# Download mouse reference (GRCm38) - 2020-A here.
# It contains the reference genome (fasta/genome.fa) and reference annotation (genes/genes.gtf).

wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"



# ============== Manage transcriptomic raw data with cellranger. ==============

# ---Install cellranger.---
wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1722198179&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=g6C6sjpkqHVB3snI-XlMFs60gseSDa8bEULI58nXbXSlb5cQ6noa-NAe0AuGuZu7CD~1OwxIhPLqSgCeZJN4asK0obM1mbwPpWWsUaLU8-AeBqV43fjQVoHvugk5IUfzsQttCTokw3Go0DGjifX0Cky-BjAUEDxuTYrr0vdSRGT7kcJ4vtnPwNCfvGOpTUPoqlXBnzJL16UcROvMen-5dyAk1TGI4H5JWUcIiI3bVfusrPkfe2kIS1NmY6iPbSZVnu2gkdCkUdU5mE26uhz2adFTGmLTc5HaftVWjOM10DARTERB7N6Vykp9~p~uHOIv5rSqr86kmNXgCxTf95MshA__"

tar -zxvf cellranger-8.0.1.tar.gz
cd ./cellranger-8.0.1/

export PATH=/home/data/SpermatogenicCells/yard/apps/cellranger-8.0.1:$PATH
which cellranger   # /home/data/SpermatogenicCells/yard/apps/cellranger-8.0.1/cellranger


# ---Alignment and gene expression quantification.---
# Map raw reads to a reference genome.
# Input: FASTQ files (.fastq, raw reads) and refdata-gex-mm10-2020-A (genes/genes.gtf).
# Run 'cellranger count'.
# Output: bam files (.bam, mapped reads), matrix files (filtered_feature_bc_matrix, filtered feature-barcode matrix), etc., that will be used later.

cellranger

# B6 cells.
cellranger count \
   --id=/home/data/SpermatogenicCells/aligned/do17815 \
   --transcriptome=/home/data/SpermatogenicCells/refdata-gex-mm10-2020-A \
   --create-bam=true \
   --fastqs=/home/data/SpermatogenicCells/rawdata/do17815/ \
   --sample=do17815 \
   --localcores=48 \
   --localmem=200

cellranger count \
   --id=/home/data/SpermatogenicCells/aligned/do17816 \
   --transcriptome=/home/data/SpermatogenicCells/refdata-gex-mm10-2020-A \
   --create-bam=true \
   --fastqs=/home/data/SpermatogenicCells/rawdata/do17816/ \
   --sample=do17816 \
   --localcores=48 \
   --localmem=200


#!!!!!!!!!!!!!!!!!!!!   
# --EXPECT_CELLS=8000
# Default is 3000. Set the number if the number of cells was preliminarily known in previous wet experiments. More accurate results.


# Cells of P5, P10, P20, P30, and P35.
**Script**: [`cellrangerCount_P5_P35.sh`](scripts/1_2_cellrangerCount_P5_P35.sh)



# ============== Manage bam files with velocyto. ==============

# ---Install velocyto and its dependencies.---
conda list velocyto

conda create -n velocyto
conda activate velocyto

conda install numpy scipy cython numba matplotlib scikit-learn h5py click
conda install pysam

conda install bioconda::velocyto.py
conda list velocyto.py
# velocyto.py               0.17.17         py311h1abe8b6_7    bioconda

# ---Spliced/unspliced/ambiguous counts quantification (or RNA velocity quantification of counts) and metadata integration.---
# Get Repeat Masker GTF file to mask out repetitive regions.
gzip -d repeat_rmsk.gtf.gz

# Input: bam files (.bam), refdata-gex-mm10-2020-A (genes.gtf), and repeat_rmsk.gtf.
# Run 'velocyto run'.
# Output: loom files (.loom) that will be used later for downstream analysis (velocyto.R/scVelo, Seurat,etc.)

velocyto --help

conda activate velocyto

velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17815/outs /home/data/SpermatogenicCells/aligned/do17815/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17816/outs /home/data/SpermatogenicCells/aligned/do17816/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do26386/outs /home/data/SpermatogenicCells/aligned/do26386/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do26387/outs /home/data/SpermatogenicCells/aligned/do26387/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17821/outs /home/data/SpermatogenicCells/aligned/do17821/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do18195/outs /home/data/SpermatogenicCells/aligned/do18195/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17824/outs /home/data/SpermatogenicCells/aligned/do17824/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do18196/outs /home/data/SpermatogenicCells/aligned/do18196/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17825/outs /home/data/SpermatogenicCells/aligned/do17825/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17827/outs /home/data/SpermatogenicCells/aligned/do17827/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
# done.


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# -b barcodes.tsv 
# -b /home/data/SpermatogenicCells/aligned/do17815/outs/filtered_feature_bc_matrix/barcodes.tsv
# According to the official docs, if you omit -b, velocyto will process every barcode it encounters. This can make the run 10–100× slower, use gigabytes of memory, and produce a lot of noise.
# Only SmartSeq2 data (one BAM per cell, usually without UMI and cell barcode tags) does not require the -b parameter.
# It is a 10X genmics single-cell RNA-Seq dataset here.

import loompy

# B6 cells.
files =  ["/home/data/SpermatogenicCells/aligned/do17815/outs/possorted_genome_bam_WJX2O.loom", "/home/data/SpermatogenicCells/aligned/do17816/outs/possorted_genome_bam_O2D4B.loom"]
loompy.combine(files,"/home/data/SpermatogenicCells/B6_merged.loom",key="Accession")

# All cells of P5, P10, P20, P30, P35, and B6.
files =  ["/home/data/SpermatogenicCells/aligned/do17815/outs/possorted_genome_bam_WJX2O.loom", "/home/data/SpermatogenicCells/aligned/do17816/outs/possorted_genome_bam_O2D4B.loom",
          "/home/data/SpermatogenicCells/aligned/do26386/outs/possorted_genome_bam_4DYG5.loom", "/home/data/SpermatogenicCells/aligned/do26387/outs/possorted_genome_bam_67K09.loom",
          "/home/data/SpermatogenicCells/aligned/do17821/outs/possorted_genome_bam_X4IXK.loom", "/home/data/SpermatogenicCells/aligned/do18195/outs/possorted_genome_bam_18F66.loom",
          "/home/data/SpermatogenicCells/aligned/do17824/outs/possorted_genome_bam_DA6KO.loom", "/home/data/SpermatogenicCells/aligned/do18196/outs/possorted_genome_bam_A70S7.loom",
          "/home/data/SpermatogenicCells/aligned/do17825/outs/possorted_genome_bam_8E32G.loom", "/home/data/SpermatogenicCells/aligned/do17827/outs/possorted_genome_bam_IKGDL.loom"]
loompy.combine(files,"/home/data/SpermatogenicCells/allten_merged.loom",key="Accession")


