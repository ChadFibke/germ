#!/software nextflow

// set up
nextflow.enable.dsl=2

params.r1 = "NULL"
params.r2 = "NULL"
params.prefix = "NULL"
params.outdir = "/outdir"

// import modules
include {
  fastqTrim;
  fastqStats} from '/germ/nextflow/modules/fastq.nf'

include {
  denovo_assembly;
  assembly_correction;
  assembly_summary
   } from '/germ/nextflow/modules/assembly.nf'

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


  fastqTrim(
    params.prefix,
    ch_read1,
    ch_read2,
    params.outdir)

  fastqStats(
    fastqTrim.out,
    params.outdir)

  denovo_assembly(
    fastqTrim.out,
    params.outdir)

  ch_assembly = fastqTrim.out
    .combine(denovo_assembly.out, by: 0)

  assembly_correction(
    ch_assembly,
    params.outdir)

  assembly_summary(
    assembly_correction.out,
    params.outdir)

}
