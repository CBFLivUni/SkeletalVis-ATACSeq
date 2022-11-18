process BAM_QC {
   
  label 'big_mem'
  publishDir "${params.outDir}/${sampleID}", mode: 'copy'

   input:
   tuple val(sampleID), path(bam)   

   output:
   path("*.txt")   

   script:
   """
   samtools index $bam
   samtools flagstat $bam > ${sampleID}_filtered_flagstat.txt
   samtools idxstats $bam > ${sampleID}_filtered_idxstats.txt
   """
 }

