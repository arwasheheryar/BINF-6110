#!/bin/bash

# Salmonella enterica genome assembly pipeline
# BINF*6110 Assignment 1 Part 2

set -euo pipefail

# Setup
PROJECT_DIR="$HOME/binf6110/binf6110/assignment01"
mkdir -p $PROJECT_DIR/{data,results/{qc,assembly,polishing,alignment,variants},figures}

cd $PROJECT_DIR
conda activate salmonella_assembly

# Step 1: Quality Control
NanoPlot \
    --fastq data/salmonella_reads.fastq.gz \
    --outdir results/qc \
    --prefix salmonella_ \
    --threads 4 \
    --N50 \
    --tsv_stats

# Step 2: Read Filtering
cd data
gunzip salmonella_reads.fastq.gz
filtlong --min_length 1000 --min_mean_q 10 \
    salmonella_reads.fastq > salmonella_filtered.fastq

# Step 3: Genome Assembly
cd $PROJECT_DIR
flye \
    --nano-hq data/salmonella_filtered.fastq \
    --out-dir results/assembly \
    --threads 8 \
    --genome-size 5m

# Step 4: Assembly Polishing
medaka_consensus \
    -i data/salmonella_filtered.fastq \
    -d results/assembly/assembly.fasta \
    -o results/polishing \
    -t 4 \
    -m r1041_e82_400bps_sup_v4.2.0

# Step 5: Reference Alignment
minimap2 \
    -ax asm5 \
    -t 8 \
    data/salmonella_reference.fna \
    results/polishing/consensus.fasta \
    | samtools view -b - \
    | samtools sort - -o results/alignment/assembly_vs_reference.sorted.bam

samtools index results/alignment/assembly_vs_reference.sorted.bam

# Step 6: Variant Calling
bcftools mpileup -Ou -f data/salmonella_reference.fna \
    results/alignment/assembly_vs_reference.sorted.bam \
    | bcftools call -mv -Oz -o results/variants/variants.vcf.gz

tabix -p vcf results/variants/variants.vcf.gz

# Step 7: Generate Analysis Files
cd results
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' \
    variants/variants.vcf.gz > variant_positions.txt

samtools idxstats alignment/assembly_vs_reference.sorted.bam > chromosome_stats.txt

bcftools query -f '%CHROM\n' variants/variants.vcf.gz | \
    sort | uniq -c | sort -rn > variants_per_chromosome.txt

# Results Summary
samtools flagstat alignment/assembly_vs_reference.sorted.bam
bcftools stats variants/variants.vcf.gz | grep "^SN"
cat assembly/assembly_info.txt
