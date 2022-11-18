process BIGWIG_MERGE {

  label 'big_mem'

  publishDir "${params.outDir}/bigwig", mode: 'copy'

   input:
   path (bigwigs)
   path (sizes)
   path (sampleTable)

   output:
   path("*.bw")   

   script:
   """

  Rscript $params.scriptDir/bigWigMerge.R \
  --bigWigs "$bigwigs"  \
  --chrSizes "$sizes"  \
  --sampleTable "$sampleTable" 
   """
 }

