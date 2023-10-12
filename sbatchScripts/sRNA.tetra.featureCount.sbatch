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


featureCounts -t exon -g Parent -O -M -a references/tetra/./L1/CombinedGff_9033_9034_WithOrthoID.all.corr.gff -o counts/tetra/L1/L1.featureCounts.tsv bams/tetra/L1/*.sam
featureCounts -t exon -g Parent -O -M -a references/tetra/./L10/CombinedGff_10707_10708_WithOrthoID.all.corr.gff -o counts/tetra/L10/L10.featureCounts.tsv bams/tetra/L10/*.sam
featureCounts -t exon -g Parent -O -M -a references/tetra/./L6/combined2508_2509.all.withOrthoID.gff -o counts/tetra/L6/L6.featureCounts.tsv bams/tetra/L6/*.sam
