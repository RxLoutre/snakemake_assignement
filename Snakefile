## --- Snakemake assignement - Code name : Nice Otter ---#
##
## Authors : @roxaneb
##

# --- Config and variables --- #

configfile: "config.json"
RESSOURCES = config['ressources']
REFERENCE = RESSOURCES + config['reference']
INDEX = RESSOURCES + config['index']
THREADS = config['threads']
RAWDATA = config['rawdata']
OUTPUT = config['output']
ANNOTATION = config['annotation']
SAMPLES = config['samples']

# --- Rules --- #

rule all:
    input:
        align_cram = expand(OUTPUT + "mapped_reads/{sample}.cram",sample = SAMPLES)
        covstats = expand(OUTPUT + "stats/{sample}_covstats.csv",sample = SAMPLES)

## help               : Prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"

## bowtie2_map        : Align pair end reads against GrCh38.
##                     /!\ The reference is expected to already be indexed for 
##                     bowtie2 and the index to be stored in the proper location
rule bowtie2_map:
  input:
    R1 = expand(RAWDATA + "{sample}_R1.fastq.gz",sample = SAMPLES),
    R2 = expand(RAWDATA + "{sample}_R2.fastq.gz",sample = SAMPLES)
  output:
    align_bam = expand(OUTPUT + "mapped_reads/{sample}.bam",sample = SAMPLES)
  params:
    index = INDEX,
    threads = THREADS
  shell:
    """
    bowtie2 -p {params.threads} -x {params.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align_bam}
    """

## bam2cram           : Copy the bam alignement into cram format
rule bam2cram:
  input:
    align_bam = expand(OUTPUT + "mapped_reads/{sample}.bam",sample = SAMPLES),
    genome_reference = REFERENCE
  output:
    align_cram = expand(OUTPUT + "mapped_reads/{sample}.cram",sample = SAMPLES)
  shell:
    """
    samtools view -C -T {input.genome_reference} -o {output.align_cram} {input.align_bam}
    """

## coverage_stat      : Produce basic coverage statistics per gene on extended exonic regions
##                      Extended = 6nt upstream and downstrem of each exon
rule coverage_stat:
  input:
    align_bam = expand(OUTPUT + "mapped_reads/{sample}.bam",sample = SAMPLES),
    genome_reference = REFERENCE,
    annotation = ANNOTATION
  output:
    covstats = expand(OUTPUT + "stats/{sample}_covstats.csv",sample = SAMPLES)
  shell:
  """

  """