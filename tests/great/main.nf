#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outDir ="data/${params.accessionNumber}"
params.scriptDir = "$baseDir/../../scripts"


include { GREAT } from '../../modules/great'

workflow test_GREAT {

  peaks = file("tests/testData/diffPeakTest.txt")

  
  GREAT(peaks,50)
}
