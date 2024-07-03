#!/bin/bash

# run using ./01.00_align-trim.sh
path="/nifH_microdiversity"
cd $path

#create output directory
mkdir $path/phylogenetic_analysis

# align with muscle 5.1
#tricho align
muscle -align $path/data_files/Tricho_top100_context.fasta -output /$path/phylogenetic_analysis/Tricho_top100_aligned.fasta
#ucyna align
muscle -align $path/data_files/UCYNA_top100_context.fasta -output /$path/phylogenetic_analysis/UCYNA_top100_aligned.fasta

# trim with trimal
#tricho trim
trimal -in $path/phylogenetic_analysis/Tricho_top100_aligned.fasta -out $path/phylogenetic_analysis/Tricho_top100_trimmed.fasta -gappyout
#ucyna trim
trimal -in $path/phylogenetic_analysis/UCYNA_top100_aligned.fasta -out $path/phylogenetic_analysis/UCYNA_top100_trimmed.fasta -gappyout
