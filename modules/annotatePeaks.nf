process ANNOTATEPEAKS {

  label 'big_mem'

  publishDir "${params.outDir}", mode: 'copy'

   input:
   path(diffPeaks)

   output:
   path("*annotated.bed")

   script:
   """
   Rscript $params.scriptDir/annotatePeaks.R \
  --diffPeaks "$diffPeaks" --annotated "${diffPeaks.baseName}_annotated.bed"
   """
 }

