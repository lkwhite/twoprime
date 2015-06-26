#! /usr/bin/env bash

#BSUB -J align
#BSUB -e log/align.%J.err
#BSUB -o log/align.%J.out
#BSUB -R "span[hosts=1]"
#BSUB -n 8

set -o nounset -o pipefail -o errexit -x

IDX=$HOME/ref/genomes/sacCer1/sacCer1
PROJECT=$HOME/projects/twoprime-seq
DATA=$PROJECT/data/2015-06-28
RESULTS=$PROJECT/results/2015-06-28
FASTQ=$DATA/SRR1158613.fastq.gz

sample=$(basename $FASTQ .fastq.gz)
bowtie2 --local -x $IDX -U $FASTQ -p 8 \
    | samtools view -F4 -bhu - \
    | samtools sort -o - $RESULTS/$sample.temp -m 8G \
    > "$sample.bam"
