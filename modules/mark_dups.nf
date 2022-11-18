process MARK_DUPS {

  label 'long_big_mem'

   input:
   tuple val(sampleID), path(bam)   

   output:
   tuple val(sampleID), path("${sampleID}_noMT_noDups.bam")   

   shell:
   '''
   samtools index !{bam}
   mkdir tmpdir
   java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=!{bam} O=!{sampleID}_noMT_noDups.bam M=!{sampleID}_dups.txt REMOVE_DUPLICATES=true TMP_DIR=tmpdir
    
   '''
 }

