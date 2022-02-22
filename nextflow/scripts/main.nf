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
  assembly_correction} from '/germ/nextflow/modules/assembly.nf'


include {
  assembly_summary;
  taxonomic_classification;
  assembly_completion} from '/germ/nextflow/modules/validation.nf'

include {
  annotation} from '/germ/nextflow/modules/annotation.nf' 

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
  // Collect PE reads for input 
  ch_read1 = Channel.fromPath( "/data/${params.r1}",
                                checkIfExists: true)

  ch_read2 = Channel.fromPath( "/data/${params.r2}",
                                checkIfExists: true)

  // Step 1: Process and summarize fastq reads
  fastqTrim(
    params.prefix,
    ch_read1,
    ch_read2,
    params.outdir)

  fastqStats(
    fastqTrim.out,
    params.outdir)

  // Step 2: Conduct assembly, correction and summary of genome  
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

  // Step 3: Check genome's taxonomy and completeness
  taxonomic_classification(
    assembly_correction.out,
    params.outdir)

  assembly_completion(
    assembly_correction.out,
    params.outdir)

  // Step 4: Assign mlst and vir/res gene presence 
  annotation(
    assembly_correction.out,
    params.outdir)

}
