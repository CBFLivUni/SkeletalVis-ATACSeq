#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "$baseDir/../../scripts"


include { BIGWIG_MERGE } from '../../modules/bigwig_merge'

workflow test_BIGWIG_MERGE {


  bigwigs = tuple(file("tests/testData/test1.bw"),
    file("tests/testData/test2.bw"))

  sizes = file("tests/testData/chrSizes.txt")

  sampleTable = file("tests/testData/bigwig_merge_sampleTable.txt")
  
 
  BIGWIG_MERGE(bigwigs,sizes,sampleTable)
}
