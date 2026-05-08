#!/bin/sh -l

#SBATCH --nodes=1                 ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 12		           ##number of requested CPUs
#SBATCH --time=04:00:00             ##time requested
#SBATCH --job-name qiime_manifest    ##name of job
#SBATCH -A microbiome  -p cpu -q standby               ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist

cd $RCAC_SCRATCH
cd buttermilk/data

pwd

#step 2, Load Modules

module load biocontainers
module load qiime2
echo "module loaded"

#step 3, Run the code you want

echo "running the code"

qiime demux subsample-paired \
  --i-sequences demux.qza \
  --p-fraction 0.3 \
  --o-subsampled-sequences demux-subsample.qza

##De-noise. Error corrections.

echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 220 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza

#Convert .qza to .qzv format

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

########FeatureTable and FeatureData summaries
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

########Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

######Alpha and beta diversity analysis (not include chao1)
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 22997 \
--m-metadata-file metadata.tsv \
--output-dir core-metrics-results



echo "End time"
date +"%d %B %Y %H:%M:%S"

#to run: sbatch qiime_pipeline.slurm
