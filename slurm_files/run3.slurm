#!/bin/bash

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH -p cpu
#SBATCH --mem=20G
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name buttermilk          ##name of job
#SBATCH -A microbiome                 ##name of the queue you are using.

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd /scratch/bell/jsenarat/buttermilk/data/run3

pwd

#step 2, Load Modules

echo "loading modules"
module load biocontainers

echo "loading qiime2"
module load qiime2

qiime feature-classifier classify-sklearn \
--i-classifier silva-classifier-v3v4.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file metadata.tsv \
--o-visualization taxa-bar-plots.qzv

qiime composition ancombc \
 --i-table table.qza \
 --m-metadata-file metadata.tsv \
 --p-formula 'Group' \
 --p-reference-levels "Group::Control diet" \
 --o-differentials ancombc-group.qza

qiime composition da-barplot \
  --i-data ancombc-group.qza \
  --p-significance-threshold 0.01 \
  --o-visualization da-barplot-group.qzv

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-l6.qza

qiime composition ancombc \
  --i-table table-l6.qza \
  --m-metadata-file metadata.tsv \
  --p-formula 'Group' \
  --p-reference-levels "Group::Control diet" \
  --o-differentials l6-ancombc-group.qza

qiime composition da-barplot \
  --i-data l6-ancombc-group.qza \
  --p-significance-threshold 0.01 \
  --p-level-delimiter ';' \
  --o-visualization l6-da-barplot-group.qzv


echo "End time"
date +"%d %B %Y %H:%M:%S"
