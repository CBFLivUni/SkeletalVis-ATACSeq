process UNIQUELY_MAPPED {

  label 'big_mem'

  publishDir "${params.outDir}/${sampleID}", mode: 'copy'

   input:
   tuple val(sampleID), path(bam)   

   output:
   tuple val(sampleID), path("${sampleID}_noMT_noDups_unique.bam"), path("${sampleID}_noMT_noDups_unique.bam.bai"),emit: filtered_bam

   script:
   """
   samtools view -b -q 10 $bam > ${sampleID}_noMT_noDups_unique.bam
   samtools index ${sampleID}_noMT_noDups_unique.bam
   """
 }

