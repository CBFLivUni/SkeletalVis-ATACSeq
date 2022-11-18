process HOMER_FIND_MOTIFS {

  label 'multi_big_mem'
  publishDir "${params.outDir}/homer", mode: 'copy'

  input:
  path(bed)

  output:
  path("SigPeaksUP")
  path("SigPeaksDOWN")
  path("targetUp.fa") , emit: targetUp
  path("backgroundUp.fa") , emit: backgroundUp
  path("targetDown.fa") , emit: targetDown
  path("backgroundDown.fa") , emit: backgroundDown

  script:
  """
   Rscript $params.scriptDir/filterPeaks.R \
  --bed "$bed" \
  --topN $params.homerTopN

  mkdir homer
  mkdir tmp
  findMotifsGenome.pl SigPeaksUP.bed hg38 SigPeaksUP -size given -dumpFasta -p $params.cpuCores -preparsedDir tmp
  mv SigPeaksUP/target.fa targetUp.fa
  mv SigPeaksUP/background.fa backgroundUp.fa

  findMotifsGenome.pl SigPeaksDOWN.bed hg38 SigPeaksDOWN -size given -dumpFasta -p $params.cpuCores -preparsedDir tmp
  mv SigPeaksDOWN/target.fa targetDown.fa
  mv SigPeaksDOWN/background.fa backgroundDown.fa

  """
}



