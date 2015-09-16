#! /usr/bin/env bash

#BSUB -J snakemaster
#BSUB -o results-twoprime/log/snakemaster_log.%J.out
#BSUB -e results-twoprime/log/snakemaster_log.%J.err

set -o nounset -o pipefail -o errexit -x

workdir="results-twoprime"
if [[ ! -d $workdir/log ]]; then
    mkdir -p $workdir/log
fi

args=' -q normal -n {threads} -o {log}.out -e {log}.err -J {params.job_name}'
snakemake --drmaa "$args" --jobs 32 --configfile ribometh-seq-config.yml
