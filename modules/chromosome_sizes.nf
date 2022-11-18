process CHROMOSOME_SIZES {

  label 'big_mem'

  publishDir "${params.outDir}", mode: 'copy'

   input:
   path(fasta)

   output:
   path("*.sizes")

   script:
   """
   gzip -df $fasta
   samtools faidx ${fasta.baseName}
   cut -f 1,2 ${fasta.baseName}.fai > ${fasta.baseName}.sizes
   """
 }