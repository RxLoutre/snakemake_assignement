rule bowtie2_map:
  input:
    index
    R1 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R1.fastq"
    R2 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R2.fastq"
  output:
    align_bam = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.bam"
  shell:
    "bowtie2 -x {input.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align_bam}"

rule bowtie2_index:
  input:
    hg19_fasta = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
  output:
    index = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/index/bowtie2/GCF_000001405.39_GRCh38.p13_genomic"
  shell:
    "bowtie2-build {input.hg19_fasta} {output.index}"

rule bam2cram:
  input:
    align_bam
  output:
    align_cram = "mapped_reads/{sample}.cram"