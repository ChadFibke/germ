#!/software nextflow

// set up
nextflow.enable.dsl=2

// Set params
params.r1 = "NULL"
params.r2 = "NULL"
params.prefix = "NULL"
params.outdir = "NULL"

ch_read1 = Channel.fromPath( "${params.r1}",
                              checkIfExists: true)

ch_read2 = Channel.fromPath( "${params.r2}",
                              checkIfExists: true)

// import modules
include {fastqStats} from '/germ/nextflow/modules/fastq.nf'


println """\
                      ~ GERM ~
            ===================================
            File Prefix      : ${params.prefix}
            Using read1      : ${params.r1}
            Using read2      : ${params.r2}
            Specified output : ${params.outdir}
            """
            .stripIndent()

workflow {
    fastqStats( params.prefix,
                ch_read1,
                ch_read2,
                params.outdir)
}
