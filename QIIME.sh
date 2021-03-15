!#/bin/bash

## Activate QIIME conda environment
conda activate qiime2-2020.11

## Create QIIME sequences artefact
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path Sequences \
   --output-path QIIME/emp-paired-end-sequences.qza
   
qiime demux emp-paired \
   --m-barcodes-file QIIME/sample-metadata.tsv \
   --m-barcodes-column barcode-sequence \
   --p-rev-comp-mapping-barcodes \
   --i-seqs QIIME/emp-paired-end-sequences.qza \
   --o-per-sample-sequences QIIME/demux-full.qza \
   --o-error-correction-details QIIME/demux-details.qza
   
qiime dada2 denoise-paired \
   --i-demultiplexed-seqs QIIME/demux-full.qza \
   --p-trim-left-f 13 \
   --p-trim-left-r 13 \
   --p-trunc-len-f 150 \
   --p-trunc-len-r 150 \
   --o-table QIIME/table.qza \
   --o-representative-sequences QIIME/rep-seqs.qza \
   --o-denoising-stats QIIME/denoising-stats.qza
   
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences QIIME/rep-seqs.qza \
  --o-alignment QIIME/aligned-rep-seqs.qza \
  --o-masked-alignment QIIME/masked-aligned-rep-seqs.qza \
  --o-tree Data/unrooted-tree.qza \
  --o-rooted-tree Data/rooted-tree.qza
  
qiime feature-classifier classify-sklearn \
   --i-classifier QIIME/silva-138-99-515-806-nb-classifier.qza \
   --i-reads QIIME/rep-seqs.qza \
   --o-classification Data/taxonomy.qza
