#!/bin/bash -l
#SBATCH -A snic2022-22-1096
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J bowtie
#SBATCH -t 24:00:00
#SBATCH --error=reports/R-%x.%j.out.StdErr
#SBATCH --output=reports/R-%x.%j.out.StdOut


# Load modules
module load bioinfo-tools
module load bowtie

# Go to correct directory
cd /proj/sllstore2017101/nobackup/b2016278/private/Neurospora/sRNA


# Run cutadapt on files


bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR4448063.cutadapt.fastq bams/crassa/SRR4448063.cutadapt.sam >bams/crassa/SRR4448063.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR751454.cutadapt.fastq bams/crassa/SRR751454.cutadapt.sam >bams/crassa/SRR751454.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR756251.cutadapt.fastq bams/crassa/SRR756251.cutadapt.sam >bams/crassa/SRR756251.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR769610.cutadapt.fastq bams/crassa/SRR769610.cutadapt.sam >bams/crassa/SRR769610.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR769748.cutadapt.fastq bams/crassa/SRR769748.cutadapt.sam >bams/crassa/SRR769748.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR4448064.cutadapt.fastq bams/crassa/SRR4448064.cutadapt.sam >bams/crassa/SRR4448064.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR755946.cutadapt.fastq bams/crassa/SRR755946.cutadapt.sam >bams/crassa/SRR755946.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/cutadapt/SRR756354.cutadapt.fastq bams/crassa/SRR756354.cutadapt.sam >bams/crassa/SRR756354.cutadapt.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR4448063.fastq bams/crassa/SRR4448063.sam >bams/crassa/SRR4448063.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR751454.fastq bams/crassa/SRR751454.sam >bams/crassa/SRR751454.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR756251.fastq bams/crassa/SRR756251.sam >bams/crassa/SRR756251.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR769610.fastq bams/crassa/SRR769610.sam >bams/crassa/SRR769610.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR769748.fastq bams/crassa/SRR769748.sam >bams/crassa/SRR769748.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR4448064.fastq bams/crassa/SRR4448064.sam >bams/crassa/SRR4448064.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR755946.fastq bams/crassa/SRR755946.sam >bams/crassa/SRR755946.summary
bowtie -q -v 0 -k 10 -S -t -p 4 references/crassa/Bowtie/ref2498_crassa reads/crassa/SRR756354.fastq bams/crassa/SRR756354.sam >bams/crassa/SRR756354.summary
