#!/bin/bash
#SBATCH --job-name=sra_getfastqfiles
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtz25002@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

#################################################################
# Download fastq files from SRA 
#################################################################

# load software
module load parallel/20180122
module load sratoolkit/3.0.1

# Dataset: GSE256438 - Purkinje cell heterogeneity in mouse cerebellum
# Study: Khouri-Farah et al. (2025) Nature Neuroscience
# FOXP genes regulate Purkinje cell diversity and cerebellar morphogenesis
# BioProject: PRJNA1079673

OUTDIR=../../data/fastq
TMPDIR=../../tmp
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR}

ACCLIST=../../metadata/subset.txt

export TMPDIR=${TMPDIR}

# use parallel to download 2 accessions at a time. 
cat $ACCLIST | parallel -j 3 "fasterq-dump -O ${OUTDIR} --split-files --temp ${TMPDIR} {}"

# compress the files 
ls ${OUTDIR}/*fastq | parallel -j 6 gzip
