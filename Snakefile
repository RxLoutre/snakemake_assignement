rule bowtie2_map:
  input:
    index
    R1 = "data/samples/A_R1.fastq"
    R2 = "data/samples/A_R2.fastq"
  output:
    align = "mapped_reads/A.bam"
  shell:
    "bowtie2 -x {input.index} -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align}"

rule bowtie2_index:
  input:
    hg19_fasta = "references/Homo_sapiens_HG19/genome/GCF_000001405.39_GRCh38.p13_genomic.fna"
  output:
    index = "references/Homo_sapiens_HG19/index/bowtie2/GCF_000001405.39_GRCh38.p13_genomic"
  shell:
    "bowtie2-build {input.hg19_fasta} {index}"