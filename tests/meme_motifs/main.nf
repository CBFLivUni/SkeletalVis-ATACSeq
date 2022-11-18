#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outDir ="data/${params.accessionNumber}"
params.motifDB = "$baseDir/../../references/db/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"
params.nmotifs = 1

include { MEME_MOTIFS } from '../../modules/meme_motifs'

workflow test_MEME_MOTIFS {

  input = file("tests/testData/target.fa")
  bg = file("tests/testData/background.fa")
   
  MEME_MOTIFS(input,bg,params.nmotifs)
}
