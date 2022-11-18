process GREAT {
  label 'big_mem'
  publishDir "${params.outDir}/GREAT", mode: 'copy'

  input:
  path(bed)
  val(topN)

  output:
  path("*.txt")

  script:
  """
 Rscript $params.scriptDir/great.R \
  --bed "$bed" \
  --topN $topN
  """
}



