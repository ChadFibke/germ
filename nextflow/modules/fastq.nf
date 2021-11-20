
// Fastq Processing

process fastqStats {
    /**
    Will tabulate the:
      - Number of reads
      - min, average, max length of reads
      - % of bases with >= Q20 scores and >= Q30
    */

    conda '/germ/envs/fastq_process.yml'

    input:
    tuple val(sample_id), file(read1), file(read2)

    output:
    stdout

    script:
    """
    seqkit stats \
        --tabular \
        --all \
        ${read1}

    seqkit stats \
        --tabular \
        --all \
        ${read2}

    """
}
