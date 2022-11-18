process DIFFBIND {

  label 'multi_big_mem'

  publishDir "${params.outDir}/DiffBind", mode: 'copy'

   input:
   path(narrowPeaks)
   path(bams)
   path (sampleTable)
   val(consensus)
   val(comparison)
   val(numerator)
   val(denominator)
  

   output:
   path ("*.png")
   path ("*PeakTable.txt"), emit: diffPeaks
   path ("*.RDS")

  

   script:
   """
   Rscript $params.scriptDir/diffBind.R \
  --peaks "$narrowPeaks" \
  --bams "$bams" \
  --sampleTable $sampleTable \
  --diffPeakTable "diffPeakTable.txt" \
  --cpuCores $params.cpuCores \
  --comparison $comparison \
  --numerator $numerator \
  --denominator $denominator \
  --useConsensus $consensus
   """
 }

