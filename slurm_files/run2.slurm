#!/bin/bash

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 4		           ##number of requested CPUs
#SBATCH -p cpu
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

#step 3, Run the code you want

echo "initiating qiime2 pipeline"

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 22997 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results


#ENTER COMMMAND TO CALCULATE CHAO1
qiime diversity alpha \
  --i-table table.qza \
  --p-metric chao1 \
  --o-alpha-diversity core-metrics-results/chao1_vector.qza


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/chao1_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/chao1-group-significance.qzv


#qiime diversity alpha-correlation \
#  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
#  --m-metadata-file metadata.tsv \
#  --o-visualization faith-pd-correlation.qzv

#qiime diversity alpha-correlation \
#  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
#  --m-metadata-file metadata.tsv \
#  --o-visualization shannon-correlation.qzv

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 40000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-Group-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column subject \
  --o-visualization core-metrics-results/bray_curtis-Group-significance.qzv \
  --p-pairwise

#For the classifier, please see run3.slurm

echo "End time"
date +"%d %B %Y %H:%M:%S"
