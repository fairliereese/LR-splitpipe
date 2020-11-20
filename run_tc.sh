#!/bin/bash
#$ -q som,bio
#$ -pe one-node-mpi 10
#$ -R y
#$ -N tc
#$ -m ea
#$ -cwd
#$ -j y

module load samtools

sam=$1
pref="${sam%.sam}"
sam_noscaff=${pref}_sorted_no_scaff.sam

# references
tc_path=/dfs6/pub/freese/mortazavi_lab/bin/TranscriptClean/
genome=/dfs6/pub/freese/mortazavi_lab/ref/mm10/mm10.fa
sjs=/data/users/freese/mortazavi_lab/ref/mm10/mm10_SJs.tsv


# remove reads that mapped to scaffold chromosomes and sort
grep -v '^@' $sam | awk ""'{if($3 !~ "_") print $0}'" "| \
          cat <(samtools view -H $sam) - | samtools view -bS | \
          samtools sort | samtools view -h > $sam_noscaff 

# run tc
time python ${tc_path}TranscriptClean.py \
    --sam $sam_noscaff \
    --genome $genome \
    --spliceJns $sjs \
    -t 8 \
    --canonOnly \
    --primaryOnly \
    --deleteTmp \
    --outprefix ${pref}    
