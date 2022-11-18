process FRAGSIZE_PLOT {

  label 'big_mem'

  publishDir "${params.outDir}/QC", mode: 'copy'

   input:
   tuple val(sampleID), path(bam), path(index)

   output:
   path("*png")
   path("*txt")

   script:
   """
   Rscript $params.scriptDir/atacQC.R \
  --bamFile "$bam"  \
  --fragmentPlot "${sampleID}_fragSize.png" \
  --fragmentTable "${sampleID}_fragTable.txt"

   """
 }

