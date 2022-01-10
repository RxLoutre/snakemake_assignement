rule bowtie2_map:
  input:
    ref = "references/GCF_000001405.39_GRCh38.p13_genomic.fna",
    R1 = "data/samples/A_R1.fastq"
    R2 = "data/samples/A_R2.fastq"
  output:
    align = "mapped_reads/A.bam"
  shell:
    "bowtie2 -1 {input.R1} -2 {input.R2} | samtools view -Sb - > {output.align}"
