process FASTQTOSAMPLEID {
   label 'big_mem'
   input:
   path sample
   path sampleTable

   output:
   tuple stdout , path(sample)

   script:
   """
   python ${params.scriptDir}/fastqToSampleID.py $sampleTable $sample.baseName
   """
}

