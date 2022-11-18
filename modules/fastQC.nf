process FASTQC {
  
  label 'big_mem'
  tag "$sampleId"
  publishDir "${params.outDir}/fastqc", mode: 'copy',
  saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(sampleId), path(read1), path(read2)

  output:
  path "*_fastqc.{zip,html}", emit: stats

  script:
  """
  fastqc -d . -t 2 -q $read1 $read2
  """

}
