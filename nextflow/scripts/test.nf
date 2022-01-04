#!/software nextflow

// set up
nextflow.enable.dsl=2

params.r1 = "NULL"
params.r2 = "NULL"
params.prefix = "NULL"
params.outdir = "/outdir"

// import modules
include {fastqTrim} from '/germ/nextflow/modules/fastq.nf'

// print options for user
println """\
                      ~ GERM ~
            ===================================
            File Prefix      : ${params.prefix}
            Using read1      : ${params.r1}
            Using read2      : ${params.r2}
            """
            .stripIndent()

// Conduct main pipeline

workflow {

  ch_read1 = Channel.fromPath( "/data/${params.r1}",
                                checkIfExists: true)

  ch_read2 = Channel.fromPath( "/data/${params.r2}",
                                checkIfExists: true)
    fastqTrim( params.prefix,
                ch_read1,
                ch_read2,
                params.outdir)
}
