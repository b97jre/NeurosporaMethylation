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


bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1aaP2-sex.cutadapt.fastq.gz bams/tetra/L1aaP2-sex.sam >bams/tetra/L1aaP2-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1aaP2.cutadapt.fastq.gz bams/tetra/L1aaP2.sam >bams/tetra/L1aaP2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R1P2-sex.cutadapt.fastq.gz bams/tetra/L1R1P2-sex.sam >bams/tetra/L1R1P2-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R1P1.cutadapt.fastq.gz bams/tetra/L1R1P1.sam >bams/tetra/L1R1P1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R2P2-sex.cutadapt.fastq.gz bams/tetra/L1R2P2-sex.sam >bams/tetra/L1R2P2-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R2P1.cutadapt.fastq.gz bams/tetra/L1R2P1.sam >bams/tetra/L1R2P1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R3P2-sex.cutadapt.fastq.gz bams/tetra/L1R3P2-sex.sam >bams/tetra/L1R3P2-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1R3P2.cutadapt.fastq.gz bams/tetra/L1R3P2.sam >bams/tetra/L1R3P2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1aP3-sex.cutadapt.fastq.gz bams/tetra/L1aP3-sex.sam >bams/tetra/L1aP3-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/./L1/L1.combinedgenome reads/tetra/cutadapt/L1aP1.cutadapt.fastq.gz bams/tetra/L1aP1.sam >bams/tetra/L1aP1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10aaP3-sex.cutadapt.fastq.gz bams/tetra/L10aaP3-sex.sam >bams/tetra/L10aaP3-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10aaP2.cutadapt.fastq.gz bams/tetra/L10aaP2.sam >bams/tetra/L10aaP2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R1P1-sex.cutadapt.fastq.gz bams/tetra/L10R1P1-sex.sam >bams/tetra/L10R1P1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R1P1.cutadapt.fastq.gz bams/tetra/L10R1P1.sam >bams/tetra/L10R1P1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R2P1-sex.cutadapt.fastq.gz bams/tetra/L10R2P1-sex.sam >bams/tetra/L10R2P1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R2P3-sex.cutadapt.fastq.gz bams/tetra/L10R2P3-sex.sam >bams/tetra/L10R2P3-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R2P2.cutadapt.fastq.gz bams/tetra/L10R2P2.sam >bams/tetra/L10R2P2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R3P3-sex.cutadapt.fastq.gz bams/tetra/L10R3P3-sex.sam >bams/tetra/L10R3P3-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10R3P3.cutadapt.fastq.gz bams/tetra/L10R3P3.sam >bams/tetra/L10R3P3.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10aP1-sex.cutadapt.fastq.gz bams/tetra/L10aP1-sex.sam >bams/tetra/L10aP1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L10/L10.combinedgenome reads/tetra/cutadapt/L10aP1.cutadapt.fastq.gz bams/tetra/L10aP1.sam >bams/tetra/L10aP1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6aaP1-sex.cutadapt.fastq.gz bams/tetra/L6aaP1-sex.sam >bams/tetra/L6aaP1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6aaP2.cutadapt.fastq.gz bams/tetra/L6aaP2.sam >bams/tetra/L6aaP2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R1P1-sex.cutadapt.fastq.gz bams/tetra/L6R1P1-sex.sam >bams/tetra/L6R1P1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R1P1.cutadapt.fastq.gz bams/tetra/L6R1P1.sam >bams/tetra/L6R1P1.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R2P1-sex.cutadapt.fastq.gz bams/tetra/L6R2P1-sex.sam >bams/tetra/L6R2P1-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R2P2.cutadapt.fastq.gz bams/tetra/L6R2P2.sam >bams/tetra/L6R2P2.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R3P2-sex.cutadapt.fastq.gz bams/tetra/L6R3P2-sex.sam >bams/tetra/L6R3P2-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6R3P3.cutadapt.fastq.gz bams/tetra/L6R3P3.sam >bams/tetra/L6R3P3.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6aP3-sex.cutadapt.fastq.gz bams/tetra/L6aP3-sex.sam >bams/tetra/L6aP3-sex.summary 2>&1
bowtie -q -v 0 -k 10 -S -t -p 4 references/tetra/L6/L6.combinedgenome reads/tetra/cutadapt/L6aP1.cutadapt.fastq.gz bams/tetra/L6aP1.sam >bams/tetra/L6aP1.summary 2>&1
