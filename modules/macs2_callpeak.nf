process MACS2_CALLPEAK {

  label 'big_mem'

  publishDir "${params.outDir}/MACS2", mode: 'copy'

   input:
   tuple val(sampleID), path(bam), path(index)

   output:
   path("${sampleID}/*.narrowPeak"), emit: narrowPeak
   tuple val(sampleID), path("${sampleID}/*_treat_pileup.bdg"), emit: bedGraph

   script:
   """
   macs2 callpeak \
  -t $bam \
  -n $sampleID \
  --outdir ${sampleID} \
  --bdg \
  -f BAMPE \
  --nomodel \
  -g 3049315783
    
   """
 }

