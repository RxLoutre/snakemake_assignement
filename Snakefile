## --- Snakemake assignement - Code name : Nice Otter ---#
##
## Authors : @roxaneb
##

# --- Config and variables --- #

configfile: "config.json"
REFERENCE = config['reference']
THREADS = config['threads']
RAWDATA = config['rawdata']
OUTPUT = config['output']
INDEX = config['index']

# --- Help Rules --- #

## help               : prints help comments for Snakefile
rule help:
    input: "Snakefile"
    shell:
        "sed -n 's/^##//p' {input}"

# --- Rules --- #

## bowtie2_map        : Align pair end reads against GrCh38.
##                     /!\ The reference is expected to already be indexed for 
##                     bowtie2 and the index to be stored in the proper location
rule bowtie2_map:
  input:
    R1 = RAWDATA + "{sample}_R1.fastq.gz",
    R2 = RAWDATA + "{sample}_R2.fastq.gz"
  output:
    align_bam = OUTPUT + "mapped_reads/{sample}.bam"
  params:
    index = INDEX,
    threads = THREADS
  shell:
    """
    bowtie2 -p {params.threads} -x {params.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align_bam}
    """
## bam2cram          : Copy the bam alignement into cram format
rule bam2cram:
  input:
    align_bam = OUTPUT + "/mapped_reads/{sample}.bam",
    genome_reference = REFERENCE
  output:
    align_cram = OUTPUT + "/mapped_reads/{sample}.cram"
  shell:
    """
    samtools view -C -T {input.genome_reference} -o {output.align_cram} {input.align_bam}
    """

  