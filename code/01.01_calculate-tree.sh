#!/bin/bash

# run using ./01.01_calculate-tree.sh
path="/nifH_microdiversity"
cd $path

#create output directories
mkdir $path/phylogenetic_analysis/tree
mkdir $path/phylogenetic_analysis/tree/tricho
mkdir $path/phylogenetic_analysis/tree/ucyna

# calculate trees with raxml
#tricho
raxmlHPC-PTHREADS-SSE3 -T 8 -f a -x 575019 -p 575019 -N 100 -m GTRGAMMA --HKY85 -O \
    -o AY768418.1 \
    -n tricho_nifH_tree.tre \
    -s $path/phylogenetic_analysis/Tricho_top100_trimmed.fasta \
    -w $path/phylogenetic_analysis/tree/tricho
#UCYNA
raxmlHPC-PTHREADS-SSE3 -T 8 -f a -x 575019 -p 575019 -N 100 -m GTRGAMMA --HKY85 -O \
    -o AY768418.1 \
    -n UCYNA_nifH_tree.tre \
    -s $path/phylogenetic_analysis/UCYNA_top100_trimmed.fasta \
    -w $path/phylogenetic_analysis/tree/ucyna