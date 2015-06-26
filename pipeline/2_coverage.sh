#! /usr/bin/env bash

#BSUB -J coverage
#BSUB -e log/coverage.%J.err
#BSUB -o log/coverage.%J.out

set -o nounset -o pipefail -o errexit -x

bamfile="SRR1158613.bam"
CHROMSIZE="$HOME/ref/genomes/sacCer1/sacCer1.chrom.sizes"

strands=(pos neg)
strand_args=('-strand +' '-strand -')

sides=(5p 3p)
side_args=('-5' '-3')
sample=SRR1158613

for strand_idx in ${!strands[@]}; do

    strand=${strands[$strand_idx]}
    strand_arg=${strand_args[$strand_idx]}

    side_files=""
    for side_idx in ${!sides[@]}; do

        side=${sides[$side_idx]}
        side_arg=${side_args[$side_idx]}

        outfile=$sample.$side.$strand.bg
        bedtools genomecov -ibam $bamfile -g $CHROMSIZE -bg \
            $strand_arg $side_arg \
            > $outfile 
        side_files="$side_files $outfile"
    done

    combined="$sample.$strand.combined.bg"
    bedtools unionbedg -i $side_files \
        | awk '{print $1, $2, $3, $4+$5}' \
        > $combined
done

