process BEDGRAPH_TO_BIGWIG {

  label 'big_mem'

  publishDir "${params.outDir}/bigwig", mode: 'copy'

   input:
   tuple val(sampleID), path(bedgraph)
   path (sizes)  

   output:
   path("*.bw"), emit: bigwig

   script:
   """
  bedGraphToBigWig $bedgraph $sizes ${sampleID}.bw
   """
 }

