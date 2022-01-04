
// Fastq Processing

process fastqTrim {
  /**
  Task:
    - Trim adaptor read-through
    - Use sliding window to remove low quality bases from 5'/3' end of read
    - Remove trimmed reads with less than 45 bases
  */

  conda '/germ/envs/fastq_process.yml'
  publishDir "$outdir/$prefix/read_summary", mode: "copy"

  input:
  val(prefix)
  file(r1)
  file(r2)
  val(outdir)

  output:
  tuple \
  val(prefix), \
  path("${prefix}_1P.fastq.gz"), \
  path("${prefix}_2P.fastq.gz"), \
  val(outdir)


  script:
  """
  trimmomatic PE \
      $r1 \
      $r2 \
      -baseout ${prefix}.fastq.gz \
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

    conda 'fastq_process'
    publishDir "$outdir/$prefix/read_summary", mode: "copy"

    input:
    val prefix
    file r1
    file r2
    val outdir

    output:
    file "*summary.tsv"

    script:
    """
    echo -e "Working with ${prefix}"
    seqkit stats \
        --tabular \
        --all \
        ${r1} > ${prefix}_1summary.tsv

    seqkit stats \
        --tabular \
        --all \
        ${r2} > ${prefix}_2summary.tsv
    """
}
