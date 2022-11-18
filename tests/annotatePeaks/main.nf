#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "$baseDir/../../scripts"


include { ANNOTATEPEAKS } from '../../modules/annotatePeaks'

workflow test_ANNOTATEPEAKS {

  bed = file("tests/testData/testBed.bed")
   
  ANNOTATEPEAKS(bed)
}
