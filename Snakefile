# Snakefile for twoprime analysis

__author__ = 'Jay Hesselberth <jay.hesselberth>'

workdir: '/vol3/home/jhessel/projects/twoprime/results'

# XXX
# configfile: 'test/samples.yaml'

SAMPLES = ('SRR1158613', 'SRR1158615')

# STRANDS and SIDES for the coverage rules
STRANDS = {'pos':'-strand +',
           'neg':'-strand -'}
SIDES  =  {'5p':'-5',
           '3p':'-3'}

SCORES = ('A','B','C')

rule all:
    input:
        'genomedata/all.genomedata'

rule fastq:
    input:
        'test/{sample}.sra'
    output:
        'raw/{sample}.fastq.gz'
    log: 
        'logs/fastq.{sample}.log'
    shell:
        'fastq-dump {input} {output}'

rule align:
    input:
        'raw/{sample}.fastq.gz'
    output:
        'alignment/{sample}/{sample}.bam'
    params:
        bwtindex = '/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1',
        id = '{sample}'
    threads: 8
    log: 
        'logs/align.{sample}.log'
    shell:
        "bowtie2 --local -x {params.bwtindex} -U {input} -p {threads} "
        "| samtools view -F4 -bhu - "
        "| samtools sort -o - alignment/{params.id}.temp -m 8G "
        "> {output} && samtools index {output}"

rule coverage:
    input:
        'alignment/{sample}/{sample}.bam'
    output:
        'coverage/{sample}/{sample}.{side}.{strand}.bg.gz'
    params:
        chromsize = '/vol3/home/jhessel/ref/genomes/sacCer1/sacCer1.chrom.sizes',
        side_arg = lambda wildcards: SIDES[wildcards.side],
        strand_arg = lambda wildcards: STRANDS[wildcards.strand]
    log: 
        'logs/coverage.{sample}.log'
    shell:
        "bedtools genomecov -ibam {input} -g {params.chromsize} -bg "
        "{params.strand_arg} {params.side_arg} "
        "| gzip -c > {output}"

def merge_coverage_input(wildcards):
    filebase = 'coverage/{sample}/{sample}.{side}.{strand}.bg.gz'
    return [filebase.format(sample = wildcards.sample,
                            strand = wildcards.strand,
                            side = side) for side in SIDES]

rule merge_coverage:
    input: merge_coverage_input
    output:
        'coverage/{sample}/{sample}.combined.{strand}.bg.gz'
    shell:
        "bedtools unionbedg -i {input} "
        "| awk '{{print $1, $2, $3, $4 + $5}}' "
        "| gzip -c > {output}"

def load_genomedata_input(wildcards):
    filebase = 'coverage/{sample}/{sample}.combined.{strand}.bg.gz'
    return [filebase.format(sample = wildcards.sample,
                            strand = wildcards.strand)]

rule load_genomedata:
    input: load_genomedata_input
    output:
       'genomedata/all.genomedata'
    log:
        'logs/genomedata.log'
    shell:
        "genomedata-load --verbose {input} {output}"

rule make_bigwigs:
    input:
        'coverage/{sample}.{strand}.combined.bg.gz'
    output:
        '{sample}.{strand}.score{score}.bw'

