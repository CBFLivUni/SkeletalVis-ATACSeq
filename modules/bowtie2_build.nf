process BOWTIE2_BUILD {
    tag "$fasta"
    label 'multi_big_mem'

    publishDir "${params.referenceDir}", mode: 'copy'

    input:
    path fasta

    output:
    tuple val("${params.index}/${fasta.baseName}"), path('bowtie2'), emit: index

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads $params.cpuCores $fasta bowtie2/${fasta.baseName}
    """
}