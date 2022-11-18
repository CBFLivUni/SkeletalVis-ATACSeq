process SAMBAMBA_MARKDUP {

  label 'multi_big_mem'

   input:
   tuple val(sampleID), path(bam)   

   output:
   tuple val(sampleID), path("${sampleID}_noMT_noDups.bam")   

   script:
   """
   mkdir tmpdir
   sambamba markdup -t $params.cpuCores --remove-duplicates --tmpdir tmpdir $bam ${sampleID}_noMT_noDups.bam
   """
 }

