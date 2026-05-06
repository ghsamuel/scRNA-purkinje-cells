#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=0-2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtz25002@uconn.edu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

module load cellranger/9.0.1

# Relative paths from project root
PROJDIR=/core/projects/GAP/GDA/gsamuel/scRNAseq
FASTQ=${PROJDIR}/data/fastq
OUT=${PROJDIR}/results/03_counts/01_cellranger
REF=/isg/shared/databases/cellranger/v2024A/refdata-gex-GRCm39-2024-A

SAMPLES=(SRR28065510 SRR28065511 SRR28065512)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

cd $FASTQ
rm -f ${SAMPLE}_S1_L001_R1_001.fastq.gz
rm -f ${SAMPLE}_S1_L001_R2_001.fastq.gz
ln -sf ${SAMPLE}_2.fastq.gz ${SAMPLE}_S1_L001_R1_001.fastq.gz
ln -sf ${SAMPLE}_3.fastq.gz ${SAMPLE}_S1_L001_R2_001.fastq.gz

mkdir -p $OUT
cd $OUT

cellranger count \
    --id=${SAMPLE} \
    --transcriptome=${REF} \
    --fastqs=${FASTQ} \
    --sample=${SAMPLE} \
    --create-bam=true \
    --localcores=16 \
    --localmem=60
