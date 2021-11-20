#!/software nextflow

// set up
nextflow.enable.dsl=2
params.NAME  = "NULL"

// import modules

include {fastqStats} from '/germ/nextflow/modules/fastq.nf'

read_ch = Channel.fromFilePairs("/germ/test/${params.NAME}_{1,2}.f*",
                                 checkIfExists: true,
                                 flat: true)


println """\
                      ~ GERM ~
            ===================================
            Prefix: ${params.NAME}
            """
            .stripIndent()




workflow {
    fastqStats(read_ch).view()

}
