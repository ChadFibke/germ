// Validation of de-novo assembly

process assembly_summary {
  /**
  Task:
    - Summarize Denovo assembly basic statistics including:
      length, # CONTIGS, N50, length, GC%
    - Estimate taxonomy
    - Estimate completness of genome
  */

  conda '/germ/envs/validation.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  path(CLEAN_ASSEMBLY)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}.assembly.summary.tsv")


  script:
  """
  seqkit stats --tabular \
    --all $CLEAN_ASSEMBLY > base.tsv

  echo "GC(%)" > gc_content

  seqkit fx2tab -C C -C G -l  $CLEAN_ASSEMBLY \
    | cut -f3,4,5,6 \
    | awk '{for(col=1;col<=NF;col++)\$col=(a[col]+=\$col)}END{print}' | awk '{print ( ((\$2 + \$3) / \$1) * 100) }' >> gc_content

   paste base.tsv gc_content > ./${PREFIX}.assembly.summary.tsv

  """
}

process taxonomic_classification {
  /**
  Task:
    - Estimate taxonamy of isolate's cleaned genome using mash
  */

  conda '/germ/envs/validation.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  path(CLEAN_ASSEMBLY)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}.taxonomy.tsv")


  script:
  """
  mash sketch $CLEAN_ASSEMBLY

  mash dist /db/refseq.genomes.k21s1000.msh ${CLEAN_ASSEMBLY}.msh > ${PREFIX}.tab
  sort -gk3  ${PREFIX}.tab | grep "GCF_" | head -n1 | cut -d_ -f1,2 > hit

  echo -e "SAMPLE\ttaxonomy" > ${PREFIX}.taxonomy.tsv

  tax=\$(grep -f hit /db/bac120_taxonomy_r202.tsv | awk -F";" '{print \$NF}')

  echo -e "${PREFIX}\t\$tax" >> ${PREFIX}.taxonomy.tsv
  """
}


process  assembly_completion {
  /**
  Task:
    - Estimate the completness of fixed assembly using busco
  */

  conda '/germ/envs/busco.yml'
  publishDir "$OUTDIR/$PREFIX/assembly", mode: "copy"

  input:
  tuple \
  val(PREFIX), \
  path(CLEAN_ASSEMBLY)
  val(OUTDIR)

  output:
  tuple \
  val(PREFIX), \
  path("${PREFIX}.busco.csv")


  script:
  """
  busco \
    --mode genome \
    --lineage_dataset /db/bacteria_odb10 \
    --in $CLEAN_ASSEMBLY \
    --out ${PREFIX}_busco \
    --offline

   egrep "%" ${PREFIX}_busco/short_summary.specific.bacteria_odb10.${PREFIX}_busco.text | tr -d "[:blank:]" > ${PREFIX}.busco.csv
  """
}
