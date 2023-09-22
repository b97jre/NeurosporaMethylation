#!/bin/bash -l
#SBATCH -A snic2022-22-1096
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J featureCounts
#SBATCH -t 2:00:00
#SBATCH --error=reports/R-%x.%j.out.StdErr
#SBATCH --output=reports/R-%x.%j.out.StdOut


# Load modules
module load bioinfo-tools
module load subread

# Go to correct directory
cd /proj/sllstore2017101/nobackup/b2016278/private/Neurospora/sRNA


# Run featureCount on files


featureCounts -t exon -g Parent -O -M -a references/crassa/2489.all.withOrthoID.gff -o counts/crassa.featureCounts.tsv bams/crassa/*.sam
