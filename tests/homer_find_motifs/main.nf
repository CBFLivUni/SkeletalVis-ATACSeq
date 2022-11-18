#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outDir ="data/${params.accessionNumber}"
params.scriptDir = "$baseDir/../../scripts"
params.cpuCores = 1
params.homerTopN = 50


include { HOMER_FIND_MOTIFS } from '../../modules/homer_find_motifs'

workflow test_HOMER_FIND_MOTIFS {

  bed = file("tests/testData/diffPeakTest.txt")
   
  HOMER_FIND_MOTIFS(bed)
}
