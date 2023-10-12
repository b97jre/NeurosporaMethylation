#!/bin/bash -l
#SBATCH -A snic2022-22-1096
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J cutadapt
#SBATCH -t 24:00:00
#SBATCH --error=reports/R-%x.%j.out.StdErr
#SBATCH --output=reports/R-%x.%j.out.StdOut


# Load modules
module load bioinfo-tools
module load cutadapt

# Go to correct directory
cd /proj/sllstore2017101/nobackup/b2016278/private/Neurospora/sRNA


# Run cutadapt on files


cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1aaP2-sex/L1aaP2-sex_S3_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1aaP2-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L1aaP2-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1aaP2/L1aaP2_S1_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1aaP2.cutadapt.fastq.gz >reads/tetra/cutadapt/L1aaP2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R1P2-sex/L1R1P2-sex_S8_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R1P2-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R1P2-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R1P1/L1R1P1_S5_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R1P1.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R1P1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R2P2-sex/L1R2P2-sex_S9_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R2P2-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R2P2-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R2P1/L1R2P1_S6_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R2P1.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R2P1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R3P2-sex/L1R3P2-sex_S10_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R3P2-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R3P2-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1R3P2/L1R3P2_S7_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1R3P2.cutadapt.fastq.gz >reads/tetra/cutadapt/L1R3P2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1aP3-sex/L1aP3-sex_S4_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1aP3-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L1aP3-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L1aP1/L1aP1_S2_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L1aP1.cutadapt.fastq.gz >reads/tetra/cutadapt/L1aP1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10aaP3-sex/L10aaP3-sex_S23_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10aaP3-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10aaP3-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10aaP2/L10aaP2_S21_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10aaP2.cutadapt.fastq.gz >reads/tetra/cutadapt/L10aaP2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R1P1-sex/L10R1P1-sex_S28_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R1P1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R1P1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R1P1/L10R1P1_S25_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R1P1.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R1P1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R2P1-sex/L10R2P1-sex_S29_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R2P1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R2P1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R2P3-sex/L10R2P3-sex_S30_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R2P3-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R2P3-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R2P2/L10R2P2_S26_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R2P2.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R2P2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R3P3-sex/L10R3P3-sex_S31_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R3P3-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R3P3-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10R3P3/L10R3P3_S27_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10R3P3.cutadapt.fastq.gz >reads/tetra/cutadapt/L10R3P3.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10aP1-sex/L10aP1-sex_S24_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10aP1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L10aP1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L10aP1/L10aP1_S22_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L10aP1.cutadapt.fastq.gz >reads/tetra/cutadapt/L10aP1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6aaP1-sex/L6aaP1-sex_S13_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6aaP1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L6aaP1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6aaP2/L6aaP2_S11_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6aaP2.cutadapt.fastq.gz >reads/tetra/cutadapt/L6aaP2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R1P1-sex/L6R1P1-sex_S18_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R1P1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R1P1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R1P1/L6R1P1_S15_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R1P1.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R1P1.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R2P1-sex/L6R2P1-sex_S19_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R2P1-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R2P1-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R2P2/L6R2P2_S16_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R2P2.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R2P2.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R3P2-sex/L6R3P2-sex_S20_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R3P2-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R3P2-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6R3P3/L6R3P3_S17_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6R3P3.cutadapt.fastq.gz >reads/tetra/cutadapt/L6R3P3.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6aP3-sex/L6aP3-sex_S14_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6aP3-sex.cutadapt.fastq.gz >reads/tetra/cutadapt/L6aP3-sex.cutadapt.summary
cutadapt -m 15 -M 35 --trimmed-only  -j 4 -a TGGAATTCTCGGGTGCCAAGG  /proj/sllstore2017101/nobackup/b2016278/private/rawdata/Nt_smallRNA/160826_D00458_0160_BC9N4VANXX/Sample_L6aP1/L6aP1_S12_L001_R1_001.fastq.gz  -o reads/tetra/cutadapt/L6aP1.cutadapt.fastq.gz >reads/tetra/cutadapt/L6aP1.cutadapt.summary
