## --- Snakemake assignement - Code name : Nice Otter ---#
##
## Authors : @roxaneb
##


## --- Config and variables --- #

REFERENCE = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/genome/GCF_000001405.39_GRCh38.p13_genomic.fna"

## --- Help Rules --- #

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
    R1 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R1.fastq.gz",
    R2 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R2.fastq.gz",
  output:
    align_bam = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.bam"
  params:
    index = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/index/bowtie2/GCF_000001405.39_GRCh38.p13_genomic"
  shell:
    """
    bowtie2 -p 8 -x {params.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align_bam}
    """
## bam2cram          : Copy the bam alignement into cram format
rule bam2cram:
  input:
    align_bam = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.bam",
    genome_reference = REFERENCE
  output:
    align_cram = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.cram"
  shell:
    """
    samtools view -C -T {input.genome_reference} -o {output.align_cram} {input.align_bam}
    """

  