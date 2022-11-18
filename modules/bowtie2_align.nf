process BOWTIE2_ALIGN {
    tag "$sampleID"
    label 'multi_big_mem'

    publishDir "${params.outDir}/bowtie2", mode: 'copy'

    input:
    tuple val(sampleID), path(read1), path(read2)
    tuple val(indexName), path(index)

    output:
    tuple val(sampleID), path('*.bam'), path('*.bai'), emit: bam
    tuple val(sampleID), path('*.txt'), emit: stats

    script:
    """
    bowtie2 --local \
    --very-sensitive-local \
    -p $params.cpuCores \
    -x bowtie2/$indexName \
    -1 $read1 \
    -2 $read2 \
    | samtools view -b | samtools sort > ${sampleID}.bam

    samtools index ${sampleID}.bam

    samtools idxstats ${sampleID}.bam > ${sampleID}_idxstats.txt

    samtools flagstat ${sampleID}.bam > ${sampleID}_flagstat.txt

    """
}



