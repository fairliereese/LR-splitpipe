#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 10
#$ -R y
#$ -N minimap
#$ -m ea
#$ -cwd
#$ -j y

fastq=$1
pref="${fastq%.fastq}"
sam=${pref}_mapped.sam
log=${pref}_minimap.log


module load minimap2

ref=/data/users/freese/mortazavi_lab/ref/mm10/mm10.fa

minimap2  \
    -t 10 \
    -ax splice:hq \
    -uf \
    --MD \
    $ref \
    $fastq \
    > $sam \
    2> $log
