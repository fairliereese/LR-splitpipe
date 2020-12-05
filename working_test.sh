#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 32
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -N test_dmux

fastq=/dfs6/pub/freese/mortazavi_lab/data/200908_splitseq/fl.fastq
python pacbio-splitpipe/demultiplex.py \
	-f $fastq \
	-o test/test \
	-t 16 \
	-i_file test/illumina_barcodes.txt \
	-rc 500