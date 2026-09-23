# Spermatogenesis
Analysis of public scRNA-Seq raw data of cells at different stages during mouse spermatogenesis.

## Environment requirements
- **Operating system**: Linux (public server)
- **python version**: 2.7.18 and 3.11.9
- **conda version**: Miniconda2-py27_4.8.3-Linux-x86_64
- **cellranger version**: cellranger-8.0.1
- **velocyto version**: 0.17.17
- **r version**: r4.3.1

## Download data

### Download scRNA-seq data
It is raw data here. Accession code **E-MTAB-6946**.  
Download all fasq files of ten time serials samples. Samples are:  
do26386 and do26387 (for postnatal days 5/P5).  
do17821 (for P10).  
do18195 (for P15).  
do17824 (for P20).  
do18196 (for P25).  
do17825 (for P30).  
do17827 (for P35).  
do17815 and do17816 (for Adult/B6).

```bash
curl https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6946/ -o sperm_html
grep ".gz" sperm_html > sperm.txt

cat sperm.txt | awk -v FS='"' '{print $8}' > sperm2.txt

awk '{print "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6946/" $0}' sperm2.txt > all_fasq_files.txt

wget -b -i all_fasq_files.txt
```
### Download reference data
Download mouse reference **(GRCm38) - 2020-A** here.  
It contains the reference genome (fasta/genome.fa) and reference annotation (genes/genes.gtf).

```bash
wget "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz"
```

## Manage transcriptomic raw data with cellranger
### Install cellranger 

```bash
wget -O cellranger-8.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1722198179&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=g6C6sjpkqHVB3snI-XlMFs60gseSDa8bEULI58nXbXSlb5cQ6noa-NAe0AuGuZu7CD~1OwxIhPLqSgCeZJN4asK0obM1mbwPpWWsUaLU8-AeBqV43fjQVoHvugk5IUfzsQttCTokw3Go0DGjifX0Cky-BjAUEDxuTYrr0vdSRGT7kcJ4vtnPwNCfvGOpTUPoqlXBnzJL16UcROvMen-5dyAk1TGI4H5JWUcIiI3bVfusrPkfe2kIS1NmY6iPbSZVnu2gkdCkUdU5mE26uhz2adFTGmLTc5HaftVWjOM10DARTERB7N6Vykp9~p~uHOIv5rSqr86kmNXgCxTf95MshA__"

tar -zxvf cellranger-8.0.1.tar.gz
cd ./cellranger-8.0.1/

export PATH=/home/data/SpermatogenicCells/yard/apps/cellranger-8.0.1:$PATH
which cellranger    # /home/data/SpermatogenicCells/yard/apps/cellranger-8.0.1/cellranger
```

### Alignment and gene expression quantification
Map raw reads to a reference genome.  
**Input**: FASTQ files (**.fastq**, raw reads) and refdata-gex-mm10-2020-A (**genes.gtf**).  
Run 'cellranger count'.  
**Output**: bam files (**.bam**, mapped reads), matrix files (**filtered_feature_bc_matrix**, filtered feature-barcode matrix), etc.

B6 cells.
```bash
cellranger

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
```

Cells of P5, P10, P20, P30, and P35.  
**Script**: [cellrangerCount_P5_P35.sh](scripts/cellrangerCount_P5_P35.sh)

## Manage bam files with velocyto
### Install velocyto and its dependencies
```bash
conda list velocyto

conda create -n velocyto
conda activate velocyto

conda install numpy scipy cython numba matplotlib scikit-learn h5py click
conda install pysam

conda install bioconda::velocyto.py
conda list velocyto.py
```

### Spliced/unspliced/ambiguous counts quantification and metadata integration
Get Repeat Masker GTF file to mask out repetitive regions.
```bash
gzip -d repeat_rmsk.gtf.gz
```

**Input**: bam files (**.bam**), refdata-gex-mm10-2020-A (**genes.gtf**), and **repeat_rmsk.gtf**.  
Run 'velocyto run'.  
**Output**: loom files (**.loom**) that will be used later for downstream analysis (velocyto.R/scVelo, Seurat, etc.)

B6 cells.
```bash
velocyto --help

conda activate velocyto

velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17815/outs /home/data/SpermatogenicCells/aligned/do17815/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200
velocyto run -m repeat_msk.gtf -o /home/data/SpermatogenicCells/aligned/do17816/outs /home/data/SpermatogenicCells/aligned/do17816/outs/possorted_genome_bam.bam /home/data/SpermatogenicCells/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 48 --samtools-memory 200

import loompy

files =  ["/home/data/SpermatogenicCells/aligned/do17815/outs/possorted_genome_bam_WJX2O.loom", "/home/data/SpermatogenicCells/aligned/do17816/outs/possorted_genome_bam_O2D4B.loom"]
loompy.combine(files,"/home/data/SpermatogenicCells/B6_merged.loom",key="Accession")
```

### Trajectory and Velocity analysis

Now, the matrix files (**filtered_feature_bc_matrix**), loom files (**.loom**), etc., are ready.  
Seurat will be used for *UMAP* and *trajectory* analysis.  
velocyto will be used for *velocity* analysis.  

<br>

**B6 Cells**  
[UMAP_B6.R](scripts/3_1_UMAP_B6.R)  
[UMAP_B6.png](<images/1_UMAP_Cell Clusters of B6 Cells (Adult).png>) | [Cldn11 Fabp3.png](<images/2_3_Feature Visualization of B6 Cells_Cldn11 Fabp3.png>) | [Dmrt1.png](<images/4_Feature Visualization of B6 Cells_Dmrt1.png>) | [Piwil1.png](<images/5_Feature Visualization of B6 Cells_Piwil1.png>) | [Tex21.png](<images/6_Feature Visualization of B6 Cells_Tex21.png>) | [Tnp1.png](<images/7_Feature Visualization of B6 Cells_Tnp1.png>)

[trajectory_B6.R](scripts/4_1_trajectory_B6.R)  
[trajectory_B6_curves.png](<images/8_Trajectory of B6 Cells_curves.png>) | [trajectory_B6_lineages.png](<images/9_Trajectory of B6 Cells_lineages.png>)

[velocity_B6.R](scripts/5_1_velocity_B6.R)  
<p align="center">
 Velocity of B6 Cells<br>
  <img src="images/10_Velocity of B6 Cells.png" alt="Velocity of B6 Cells" width="600">
</p>

<br>

**Cells of P5, P10, P20, P30, P35 and B6.**  
[UMAP_all.R](scripts/3_2_UMAP_all.R)  
[UMAP_all.png](<images/11_Cell Clusters of All Sample Cells (P5, P10, ..., P35, B6).png>) | [Cldn11.png](<images/12_Feature Visualization of All Sample Cells_Cldn11.png>) | [Fabp3.png](<images/13_Feature Visualization of All Sample Cells_Fabp3.png>) | [Dmrt1.png](<images/14_Feature Visualization of All Sample Cells_Dmrt1.png>) | [Piwil1.png](<images/15_Feature Visualization of All Sample Cells_Piwil1.png>) | [Tex21.png](<images/16_Feature Visualization of All Sample Cells_Tex21.png>) | [Tnp1.png](<images/17_Feature Visualization of All Sample Cells_Tnp1.png>)

[trajectory_all.R](scripts/4_2_trajectory_all.R)   
[trajectory_all_curves_noRoot.npg](<images/18_Trajectory of All Sample Cells_curves_Legend.png>) | [trajectory_all_lineages_noRoot.npg](<images/19_Trajectory of All Sample Cells_lineages_Legend.png>)  
[trajectory_all_curves_SpermatogoniaRoot.npg](<images/20_Trajectory of All Sample Cells_curves_Legend.png>) | [trajectory_all_lineages_SpermatogoniaRoot.npg](<images/21_Trajectory of All Sample Cells_lineages_Legend.png>)

[velocity_all.R](scripts/5_2_velocity_all.R)  
<p align="center">
 Velocity of All Sample Cells<br>
  <img src="images/22_Velocity of All Sample Cells.png" alt="Velocity of All Sample Cells" width="600">
</p>
