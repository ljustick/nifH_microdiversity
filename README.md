# nifH_microdiversity

[![DOI](https://zenodo.org/badge/656212449.svg)](https://zenodo.org/doi/10.5281/zenodo.12636629)


Code to accompany: Global Phylogeography and Microdiversity of the Marine Diazotrophic Cyanobacteria Trichodesmium and UCYN-A

Angie Nguyen, Lucas J. Ustick, Alyse A. Larkin, Adam C. Martiny

## QUIMEII processing (Written by Nguyen)
1) HPC_Commands_nifH_P18_AMT.sh
- This script is the QIIME2 pipeline used to trim primers from and quality filter sequence data
- Input files: 
	- raw sequence data as fastq files
	- manifest reformat: this is a file with the following information in this specific format: sample-id,filepath,direction for each fastq file
	- metadata 
- Output files:
	- dna-sequences.fasta: trimmed and quality filtered sequences with OTU id, will appear in the output directory "nifH-P18-AMT-cutadapt-rep-seqs-output"
	- PAMT-feature-table.tsv: table of features (each unique sequence with an assigned OTU id) and its abundance in each sample

2) HPC_transeq.sh
- Transeq reads nucleotide sequences in up to six different frames and writes the corresponding amino acid sequence translations. All 6 frames were run, and the most appropriate frame, yielding the most complete protein sequences, was kept for further analysis.
- Input file: dna-sequences.fasta - trimmed and quality filtered sequences with OTU id
- Output file: translated-dna-sequences - translated protein sequences with OTU id

3) HPC_PAMT_blast.sh
- Used Ensembl BLAST (Basic Local Alignment Search Tool) to match protein sequences to a known database. The goal is to assign each OTU to a particular organism.
- We curated our own database to align our query sequences against.
- Input file: translated-dna-sequences - translated protein sequences with OTU id
- Output File: nifH-PAMT-blast-output.txt - OTUs with matched organisms and relevant statistics

## Data analysis (Written by Ustick)

### Phylogenetic Tree Analysis
- **01.00_align-trim.sh** Align fasta sequences, and trim
    - align reads with muscle 5.1
    - trim sequences with trimAL 1.4.1
- Determine best substitution model using megaX
    - results can be found in phylogenetic_analysis/model***
- **01.01_calculate-tree.sh** calculate phylogenetic trees
    - calculate tree with raxml-8.2.12
- **Visualize using ITOL**
    - use the folling files
        - RAxML_bipartitionsBranchLabels.UCYNA_nifH_tree.tre
        - RAxML_bipartitionsBranchLabels.tricho_nifH_tree.tre
    - Save ITOL outputs for R analysis
        - itol_ucyna_newick.txt
        - itol_tricho_newick.txt

### Environmental Context
- **02.00_phylogeography.R** Code to analyze the clades within UCYN-A/Tricho.
    - Identify clades
    - Calculate relative abundances
    - Figure 1: sample map
    - Figure 2: top ASV's
    - Figure 3: biogeography
    - Figure 4: statistical analysis
