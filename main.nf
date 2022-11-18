nextflow.enable.dsl=2

//define the parameters and the default values

//outdirectory for the results
params.outDir = "$baseDir/data/${params.accessionNumber}"

//directory with the scripts
params.scriptDir = "$baseDir/scripts"

//sample table for the experiment
params.sampleTable = "$baseDir/params/sampleTables/${params.accessionNumber}_sampleTable.txt"

//define the number of cpus to use coverting sra to fastq
params.cpuCores=1

//directory with the fastqFiles
params.fastqFileDir = "$params.outDir/fastqFiles/*.fastq.gz"
downloadFiles = file(params.fastqFileDir).isEmpty()

//directory with the index files
params.referenceDir = "$baseDir/references"

params.fasta = "ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

//species to use for downloading the index file and enrichment analysis
params.species="human"

//index for bowtie
params.index = "$params.referenceDir/bowtie2"
params.indexName = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
makeIndex = file(params.index).isEmpty()

params.replicates = true

//directory with the motifs for meme-tomtom
params.motifDB = "$params.referenceDir/db/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme"

//should donor/tissue level consensus peaks (specifed by the DB_TISSUE column in the samplesheet)
//be used in the DiffBind analysis?
params.consensus = "FALSE"

//the maximum number of TF motifs for meme to identify
params.nmotifs = 15

//the number of top up and down altered genes for homer motif enrichment
params.homerTopN = 250

def helpMessage() {
  log.info"""
  =========================================
  Skeletalvis-ATAC Pipeline
  =========================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run ATAC.nf -profile slurm -params-file params/<accessionNumber>.yaml -with-singularity atac.img

  Required arguments:
  -profile                      Configuration profile to use. <local, slurm>
  -params-file                  Yaml file containing parameters for the analysis
  --with-singularity            Recommended to run the pipeline within the provided singularity container

  """.stripIndent()
}

params.help = false
if (params.help){
  helpMessage()
  exit 0
}

//load all the modules

include { GET_FASTQFILES_ENA } from './modules/get_fastqfiles_ena'
include { FASTQTOSAMPLEID } from './modules/fastqToSampleID'
include { FASTQC } from './modules/fastQC'
include { BOWTIE2_BUILD } from './modules/bowtie2_build'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align'
include { REMOVE_MT } from './modules/remove_mt'
include { MARK_DUPS } from './modules/mark_dups'
include { SAMBAMBA_MARKDUP } from './modules/sambamba_markdup'
include { UNIQUELY_MAPPED } from './modules/uniquely_mapped'
include { BAM_QC } from './modules/bam_qc'
include { MACS2_CALLPEAK } from './modules/macs2_callpeak'
include { CHROMOSOME_SIZES } from './modules/chromosome_sizes'
include { BEDGRAPH_TO_BIGWIG } from './modules/bedgraph_to_bigwig'
include { DIFFBIND } from './modules/diffbind'
include { FRAGSIZE_PLOT } from './modules/fragsize_plot'
include { BIGWIG_MERGE } from './modules/bigwig_merge'
include { ANNOTATEPEAKS } from './modules/annotatePeaks'
include { HOMER_FIND_MOTIFS } from './modules/homer_find_motifs'
include { MEME_MOTIFS as MEME_MOTIFS_UP  } from './modules/meme_motifs'
include { MEME_MOTIFS as MEME_MOTIFS_DOWN} from './modules/meme_motifs'

include { GREAT } from './modules/great'


//subworkflow to download the fastq files from ENA/SRA if needed
workflow DOWNLOAD_DATA {

  take:
  sampleTable

  main:
  fastqFiles = GET_FASTQFILES_ENA(sampleTable)

  emit:
  fastqFiles


}


//subworkflow to match fastq files to the sampleID and align using bowtie2
workflow ALIGN {

  take:
  sampleTable
  rawData


  main:
  

    //create the bowtie2 index if not pesent in the ref dir
    if(makeIndex){

      ref = BOWTIE2_BUILD(params.fasta)
      index = ref.index
      
      } else {

        index = tuple("bowtie2/$params.indexName",file(params.index))

      }

      //group the data if paired reads
      groupedReads = rawData.map { file ->
        def key = file.name.toString().split("_R\\d*")[0]
        return tuple(key, file)
      }
      .groupTuple().view()

      //map the fastq files to the sampleID and group
      samples_grouped = FASTQTOSAMPLEID(rawData,sampleTable)
      .groupTuple(sort:true)
      .map { id, fastqs -> tuple( id, fastqs[0], fastqs[1] ) }
      .view()

      FASTQC(samples_grouped)

      //align the reads and fitler the bams
      BOWTIE2_ALIGN(samples_grouped,index)


      emit: BOWTIE2_ALIGN.out.bam

      
    }

    //subworkflow to filter bams by mitochondria, remove duplicate, multimapped and do QC
    workflow FILTERBAM {

      take: bam

      main:

      REMOVE_MT(bam) | SAMBAMBA_MARKDUP | UNIQUELY_MAPPED | (BAM_QC & FRAGSIZE_PLOT)


      emit: UNIQUELY_MAPPED.out.filtered_bam

    }

    //subworkflow to call peaks, merge bigwigs files per condition and annotate differential peaks
    workflow FINDPEAKS {

      take: 
      bam_with_index
      sampleTable

      main:

      MACS2_CALLPEAK(bam_with_index)
      sizes = CHROMOSOME_SIZES(params.fasta)
      BEDGRAPH_TO_BIGWIG(MACS2_CALLPEAK.out.bedGraph,sizes)

      if(params.replicates) {

        bigwigs = BEDGRAPH_TO_BIGWIG.out.bigwig.collect()
        BIGWIG_MERGE(bigwigs,sizes,sampleTable)

        peaks = MACS2_CALLPEAK.out.narrowPeak | collect

        bams = bam_with_index.map { id, bam, index -> tuple( bam, index ) }.collect().view()

        DIFFBIND(peaks,bams,sampleTable,params.consensus,params.comparison,params.numerator,params.denominator)

        ANNOTATEPEAKS(DIFFBIND.out.diffPeaks)

        HOMER_FIND_MOTIFS(DIFFBIND.out.diffPeaks)
        MEME_MOTIFS_UP(HOMER_FIND_MOTIFS.out.targetUp,HOMER_FIND_MOTIFS.out.backgroundUp, params.nmotifs)
        MEME_MOTIFS_DOWN(HOMER_FIND_MOTIFS.out.targetDown,HOMER_FIND_MOTIFS.out.backgroundDown, params.nmotifs)
        GREAT(DIFFBIND.out.diffPeaks, 5000)
      }



    }



  workflow {

  //read in the sample table for the experiment
  sampleTable = file( params.sampleTable )

  if(downloadFiles){

    rawData = DOWNLOAD_DATA(sampleTable).fastqFiles.flatten()

    } else {

      rawData = channel.fromPath( params.fastqFileDir).view()

    }

  //align, filter and find the peaks
  filtered_bams = ALIGN(sampleTable,rawData) | FILTERBAM
  FINDPEAKS(filtered_bams,sampleTable)

}



