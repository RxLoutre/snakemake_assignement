genome_reference = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/genome/GCF_000001405.39_GRCh38.p13_genomic.fna"

rule bowtie2_index:
  input:
    hg19_fasta = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
  output:
    index = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/index/bowtie2/GCF_000001405.39_GRCh38.p13_genomic"
  shell:
    "bowtie2-build --threads 8 {input.hg19_fasta} {output.index}"

rule bowtie2_map:
  input:
    R1 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R1.fastq.gz",
    R2 = "/Users/roxaneboyer/Bioinformatic/data/vUMC/rawdata/{sample}_R2.fastq.gz",
    index = "/Users/roxaneboyer/Bioinformatic/ressources/references-database/Homo_sapiens_HG19/index/bowtie2/GCF_000001405.39_GRCh38.p13_genomic"
  output:
    align_bam = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.bam"
  shell:
    "bowtie2 -p 8 -x {input.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align_bam}"

rule bam2cram:
  input:
    align_bam = "/Users/roxaneboyer/Bioinformatic/data/vUMC/mapped_reads/{sample}.bam"
  output:
    align_cram = "mapped_reads/{sample}.cram"