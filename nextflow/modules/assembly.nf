// De-novo assembly and processing

process denovo_assembly {
  /**
  Task:
    - Perform Denovo assembly using trimmed reads
  */

  conda '/germ/envs/assembly.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  file(R1_TRIMMED), \
  file(R2_TRIMMED)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}.skesa.fa")


  script:
  """
  skesa \
    --cores 4 \
    --vector_percent 1 \
    --fastq $R1_TRIMMED,$R2_TRIMMED > ${PREFIX}.skesa.fa
  """
}


process assembly_correction {
  /**
  Task:
    - Align trimmed reads to Denovo assembly
    - Use alignment to correct assembled contigs
  */

  conda '/germ/envs/assembly.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  path(R1_TRIMMED), \
  path(R2_TRIMMED), \
  path(ASSEMBLY)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}.fix.fasta")


  script:
  """
  bwa index $ASSEMBLY
  samtools faidx $ASSEMBLY

  bwa mem $ASSEMBLY $R1_TRIMMED $R2_TRIMMED \
    | samtools sort -o ${PREFIX}.bam -

  samtools index ${PREFIX}.bam

  pilon \
    --genome $ASSEMBLY \
    --frags ${PREFIX}.bam \
    --fix bases \
    --mindepth 1 \
    --minmq 20 \
    --minqual 20 \
    --output ${PREFIX}.pilon

  # find low coverage contigs
  seqkit seq --name ${PREFIX}.pilon.fasta | awk -F_ '{if(\$3<10) print \$0}' > low_cov_index

  seqkit seq --min-len 499 ${PREFIX}.pilon.fasta \
  | seqkit grep --invert-match --pattern-file low_cov_index > ${PREFIX}.fix.fasta

  """
}
