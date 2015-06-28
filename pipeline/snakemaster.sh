#! /usr/bin/env bash

#BSUB -J snakemaster
#BSUB -o snakemaster.%J.out
#BSUB -e snakemaster.%J.err

args=' -q normal -n {threads} -o {log}.out -e {log}.err -J {params.job_name}'
snakemake --drmaa "$args" --jobs 32
