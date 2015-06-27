# Snakefile for twoprime analysis

__author__ = 'Jay Hesselberth <jay.hesselberth>'
__version__ = '0.01b'

SAMPLES = ('SRR1158613').split()

# STRANDS and SIDES for the coverage rules
STRANDS = {'pos':'-strand +',
           'neg':'-strand -'}
SIDES  =  {'5p':'-5',
           '3p':'-5'}

rule all:
    input:
        expand("{sample}.{strand}.signal.bw", sample=SAMPLES)
        expand("{sample}.{strand}.scoreA.bw", sample=SAMPLES)
        expand("{sample}.{strand}.scoreB.bw", sample=SAMPLES)
        expand("{sample}.{strand}.scoreC.bw", sample=SAMPLES)
        'hub/hub.txt'

rule align:
    input:
        'raw/{sample}.fastq.gz'
    output:
        'alignment/{sample}/{sample}.bam'
    params:
        job_name='align.{sample}'
        sample='{sample}'
    threads: 8
    log: 
        'logs/align.{sample}.log'
    shell:
        'bowtie2 --local -x {BWTINDEX} -U {input} -p {threads} \
            | samtools view -F4 -bhu - \
            | samtools sort -o - {sample}.temp -m 8G \
            > {output} && samtools index {output}'

rule coverage:
    input:
        bam='alignment/{sample}/{sample}.bam'
        chromsize='/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1.chrom.sizes'
        strand=expand(STRANDS)
        side=expand(SIDES)
    output:
        'coverage/{sample}/{sample}.{input.side}.{input.strand}.bg.gz'
    shell:
        'bedtools genomecov -ibam {input.bam} -g {input.chromsize} -bg \
            {strand_arg} {side_arg} \
            | gzip -c > {output}'

rule merge_coverage:
    input:
    output:
          
