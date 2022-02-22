process  annotation {
  /**
  Task:
    - Identify mlst for isolates using auto detect
    - Use abricate to identify virulence/antimicroial resistance genes from ncbi
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
  path("${PREFIX}.mlst.tsv"), \
  path("${PREFIX}.abricate.tsv")


  script:
  """
  mlst $CLEAN_ASSEMBLY > ${PREFIX}.mlst.tsv

  abricate $CLEAN_ASSEMBLY > ${PREFIX}.abricate.tsv
  """
}
