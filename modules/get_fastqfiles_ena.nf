process GET_FASTQFILES_ENA {
  label 'big_mem'
  publishDir path: "${params.outdir}/fastqFiles", mode: 'copy'


  input:
  path (sampleTable)

  output:
  path ('*.fastq.gz'), emit: fastqFiles


  script:

  """
  Rscript ${params.scriptDir}/getFastqFilesFromENA.R --accessionNumber "${params.ENA}" --sampleTable $sampleTable --downloadedFiles "downloadedFile.txt"
  """

}


