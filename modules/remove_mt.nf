process REMOVE_MT {

  label 'big_mem'

   input:
   tuple val(sampleID), path(bam), path(bamIndex)   

   output:
   tuple val(sampleID), path("${sampleID}_noMT.bam")   

   script:
   """
   samtools idxstats $bam | cut -f 1 | grep -v MT | xargs samtools view -b $bam > ${sampleID}_noMT.bam
   """
 }

