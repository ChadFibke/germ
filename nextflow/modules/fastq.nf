
// Fastq Processing

process fastqTrim {
  /**
  Task:
    - Trim adaptor read-through
    - Use sliding window to remove low quality bases from 5'/3' end of read
    - Remove trimmed reads with less than 45 bases
  */

  conda '/germ/envs/fastq_process.yml'
  publishDir "$OUTDIR/$PREFIX/read_summary", mode: "copy"

  input:
  val(PREFIX)
  file(R1)
  file(R2)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}_1P.fastq.gz"), \
  path("${PREFIX}_2P.fastq.gz")


  script:
  """
  trimmomatic PE \
      $R1 \
      $R2 \
      -baseout ${PREFIX}.fastq.gz \
      ILLUMINACLIP:/germ/dependencies/trimmomatic/Illumina_adapters.fa:1:45:15 \
      LEADING:5 \
      TRAILING:5 \
      SLIDINGWINDOW:10:20 \
      MINLEN:45 \
      TOPHRED33
  """

}


process fastqStats {
    /**
    Task:
      - Tabulate number of reads
      - Tabulate min, average, max length of reads
      - Tabulate % of bases with >= Q20 scores and >= Q30
      - Provide surface level fastqc metrics
    */

    conda '/germ/envs/fastq_process.yml'
    publishDir "$OUTDIR/$PREFIX/read_summary", mode: "copy"

    input:
    tuple \
    val(PREFIX), \
    file(R1_TRIMMED), \
    file(R2_TRIMMED)
    val(OUTDIR)

    output:
    tuple \
    val(PREFIX), \
    path("${PREFIX}.seqkit1.tsv"), \
    path("${PREFIX}.seqkit2.tsv"), \
    path("${PREFIX}.fastqc1.summary.txt"), \
    path("${PREFIX}.fastqc2.summary.txt")

    script:
    """
    echo -e "Working with ${PREFIX}"
    seqkit stats \
      --tabular \
      --all \
      $R1_TRIMMED > ${PREFIX}.seqkit1.tsv

    seqkit stats \
      --tabular \
      --all \
      $R2_TRIMMED > ${PREFIX}.seqkit2.tsv

    fastqc \
      $R1_TRIMMED \
      $R2_TRIMMED

    unzip -p ${PREFIX}_1P_fastqc.zip ${PREFIX}_1P_fastqc/summary.txt > ${PREFIX}.fastqc1.summary.txt
    unzip -p ${PREFIX}_2P_fastqc.zip ${PREFIX}_2P_fastqc/summary.txt > ${PREFIX}.fastqc2.summary.txt
    """
}
