#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "$baseDir/../../scripts"


include { FRAGSIZE_PLOT } from '../../modules/fragsize_plot'

workflow test_FRAGSIZE_PLOT {


  bam = tuple("testBam",
    file("tests/testData/testChr6_r2.bam"),
    file("tests/testData/testChr6_r2.bam.bai"))
  
 
  FRAGSIZE_PLOT(bam)
}
