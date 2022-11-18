process MEME_MOTIFS {
  label 'big_mem'
  publishDir "${params.outDir}/meme", mode: 'copy'

  input:
  path(input)
  path(background)
  val(nmotifs)

  output:
  path("$input.baseName")
  path("tomtom_${input.baseName}")

  script:
  """
  export TMPDIR="tmp/"
  
  #meme using homer background sequences
  meme $input -nmotifs $nmotifs -maxw 20 -dna -oc ${input.baseName} -objfun de -neg $background

  #tom tom to match de novo motifs to known motifs
  tomtom -no-ssc -oc tomtom_${input.baseName} -min-overlap 5 -dist pearson ${input.baseName}/meme.txt $params.motifDB

  """
}



