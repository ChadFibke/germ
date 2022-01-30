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


process assembly_summary {
  /**
  Task:
    - Summarize Denovo assembly basic statistics including:
      length, # CONTIGS, N50, length, GC%
    - Estimate taxonomy
    - Estimate completness of genome
  */

  conda '/germ/envs/annotation.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  path(CLEAN_ASSEMBLY)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}_assembly_summary.tsv")


  script:
  """
  seqkit stats --tabular \
    --all $CLEAN_ASSEMBLY > base.tsv

  echo "GC(%)" > gc_content

  seqkit fx2tab -C C -C G -l test.fix.fasta \
    | cut -f3,4,5,6 \
    | awk '{for(col=1;col<=NF;col++)\$col=(a[col]+=\$col)}END{print}' | awk '{print ( ((\$2 + \$3) / \$1) * 100) }' >> gc_content

   paste base.tsv gc_content > ./${PREFIX}_assembly_summary.tsv

  """
}
