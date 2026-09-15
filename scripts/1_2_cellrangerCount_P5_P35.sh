#!/bin/bash
# ==================================================================
# Script name: 1_2_cellrangerCount_P5_P35.sh
# Purpose: run 'cellranger count' for multiple samples
# Usage: sbatch scripts/1_2_cellrangerCount_P5_P35.sh
# Author: Xiaoya
# AI assistance: DeepSeek (for code debugging)
# Date: 2026-09-10

# Remember to make sure:
# mkdir -p /home/data/SpermatogenicCells/logs
# ls -d /home/data/SpermatogenicCells/rawdata/do26386/
# ls /home/data/SpermatogenicCells/rawdata/do26386/*.fastq.gz | head
# ==================================================================

#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=32
#SBATCH --mem=140G
#SBATCH --time=24:00:00
#SBATCH --array=1-8%2
#SBATCH --output=/home/data/SpermatogenicCells/logs/%A_%a.out
#SBATCH --error=/home/data/SpermatogenicCells/logs/%A_%a.err

set -euo pipefail

export PATH=/home/data/SpermatogenicCells/yard/apps/cellranger-8.0.1:$PATH

TRANSCRIPTOME=/home/data/SpermatogenicCells/refdata-gex-mm10-2020-A
RAW_DATA=/home/data/SpermatogenicCells/rawdata
OUT_BASE=/home/data/SpermatogenicCells/aligned

SAMPLES=("do26386" "do26387" "do17821" "do18195" \
         "do17824" "do18196" "do17825" "do17827")
SAMPLE="${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"

echo "[$(date)] START $SAMPLE"

mkdir -p "$OUT_BASE"
cd "$OUT_BASE" || { echo "ERROR: cannot cd to $OUT_BASE"; exit 1; }

if [ -d "$SAMPLE" ]; then
    echo "Output dir $SAMPLE already exists, skipping"
    exit 0
fi

trap 'rm -rf "${SAMPLE}/SC_RNA_COUNTER_CS" "${SAMPLE}/"_*' EXIT

cellranger count \
    --id="$SAMPLE" \
    --transcriptome="$TRANSCRIPTOME" \
    --create-bam=true \
    --fastqs="$RAW_DATA/$SAMPLE" \
    --sample="$SAMPLE" \
    --localcores=24 \
    --localmem=120

echo "[$(date)] DONE $SAMPLE"