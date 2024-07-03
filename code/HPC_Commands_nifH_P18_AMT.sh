#!/bin/bash
#$ -N qiime-nifH-PAMT-dadaON
#$ -o qiime-nifH-PAMT-dadaON.out
#$ -e qiime-nifH-PAMT-dada.err
#$ -m beas
#$ -pe openmp 24-64
#$ -q bio,abio,abio128
#$ -ckpt restart

module load anaconda/3.6-4.3.1
module load qiime2
source activate qiime2-2019.7

cd /pub/angien4/P18-AMT/

## Import Data
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /pub/angien4/P18-AMT/nifH_P18-AMT_manifest_reformat.csv --output-path nifH-P18-AMT-paired-end-sequences.qza --input-format PairedEndFastqManifestPhred33


## Visualize sequence quality
qiime demux summarize --i-data nifH-P18-AMT-paired-end-sequences.qza --o-visualization nifH-P18-AMT-demux.qzv

# Cutadapt: Removing Primers from sequences before denoising steps #
qiime cutadapt trim-paired --i-demultiplexed-sequences nifH-P18-AMT-paired-end-sequences.qza --p-front-f TGYGAYCCNAARGCNGA --p-front-r ADNGCCATCATYTCNCC --p-error-rate 0.2 --p-discard-untrimmed --o-trimmed-sequences nifH-P18-AMT-cutadapt-trim-seqs.qza --verbose --p-no-indels

# Visualize sequences without primers #
qiime demux summarize --i-data nifH-P18-AMT-cutadapt-trim-seqs.qza --o-visualization nifH-P18-AMT-cutadapt-trim-seqs.qzv   
  
lsQuality Filter your Sequences aka Denoising Paired-End

qiime dada2 denoise-paired --i-demultiplexed-seqs nifH-P18-AMT-cutadapt-trim-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 278 --p-trunc-len-r 243 --p-max-ee-f 4 --p-max-ee-r 4 --p-n-reads-learn 800000 --o-table nifH-P18-AMT-cutadapt-table.qza --o-representative-sequences nifH-P18-AMT-cutadapt-rep-seqs.qza --o-denoising-stats nifH-P18-AMT-cutadapt-denoising-stats.qza --p-n-threads $CORES


## Summarize Denoising/Quality Filtering Step 
qiime feature-table summarize --i-table nifH-P18-AMT-cutadapt-table.qza --o-visualization nifH-P18-AMT-cutadapt-table.qzv --m-sample-metadata-file /data.users/angien4/P18-AMT/mR143-L1/P18-AMT-metadata_new-reformat.txt 

qiime feature-table tabulate-seqs --i-data nifH-P18-AMT-cutadapt-rep-seqs.qza --o-visualization nifH-P18-AMT-cutadapt-rep-seqs.qzv

## Visualize Denoising Stats 
qiime metadata tabulate --m-input-file nifH-P18-AMT-cutadapt-denoising-stats.qza --o-visualization nifH-P18-AMT-cutadapt-denoising-stats.qzv

## Export Rep-Seqs to FastQ file for BLAST
qiime tools export --input-path nifH-P18-AMT-cutadapt-rep-seqs.qza --output-path nifH-P18-AMT-cutadapt-rep-seqs-output 

#Export table.qza into .biom and OTU table
qiime tools export --input-path /pub/angien4/P18_AMT/nifH-P18-AMT-cutadapt-table.qza --output-path /pub/angien4/P18_AMT/blast
cd blast
biom convert -i feature-table.biom -o PAMT-feature-table.tsv --to-tsv 
biom head -i PAMT-feature-table.tsv


conda deactivate
