#!/bin/bash
#$ -N nifH-C-13-transeqON
#$ -o nifH-C-13-transeqON.out
#$ -e nifH-C-13-transeq.err
#$ -m beas
#$ -pe openmp 24-64
#$ -q bio,abio,abio128
#$ -ckpt restart

module load emboss/6.5.7

cd /dfs3/nifty/processed/C13.5/qiime2/cutadapt/nifH-C13-5-cutadapt-rep-seqs-output/trans-AN

transeq -sequence /dfs3/nifty/processed/C13.5/qiime2/cutadapt/nifH-C13-5-cutadapt-rep-seqs-output/dna-sequences.fasta -outseq translated-dna-sequences -frame 2

#frame will be different because these data used different primers, run frame 6 to see all possible frames, look at sequences and then choose the best one
#you need to translate the sequences into different frames and look at the sequences, * denote stop codons

conda deactivate
