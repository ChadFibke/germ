
// Fastq Processing

process fastqStats {
    /**
    Will tabulate the:
      - Number of reads
      - min, average, max length of reads
      - % of bases with >= Q20 scores and >= Q30
    */

    conda '/germ/envs/fastq_process.yml'
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
