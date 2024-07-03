#!/bin/bash
#$ -N nifH-P18-AMT-blastpON
#$ -o nifH-P18-AMT-blastpON.out
#$ -e nifH-P18-AMT-blastp.err
#$ -m beas
#$ -pe openmp 24-64
#$ -q bio,abio,abio128
#$ -ckpt restart

module load anaconda/3.7-5.3.0

cd /pub/angien4/P18_AMT/nifH-P18-AMT-cutadapt-rep-seqs-output

#makeblastdb -in framebot.fa -out nifH-database-blastp -dbtype prot -parse_seqids
#use new database

blastp -db zehr_delm_degap_AA_db.blastp -num_threads 20 -query translated-dna-sequences -evalue 0.00001 -out nifH-PAMT-blast-output.txt -max_target_seqs 1 -outfmt 6

conda deactivate