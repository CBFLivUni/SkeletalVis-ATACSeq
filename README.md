
# Nextflow pipeline for ATAC-seq data analysis

## Introduction

This bioinformatics pipeline enables reproducible analyses of ATAC-sequencing data.

The **SkeletalVis-ATAC-seq** pipeline is built using [Nextflow](https://www.nextflow.io), a portable workflow tool to run tasks across multiple compute infrastructures. This pipeline uses singularity containes containing all the software needed to run the analysis, making installation simple and the results reproducible.

## Pipeline summary

The **SkeletalVis-ATAC-seq** pipeline takes a sample table and a parameter file defining the experiment as input. If not provided fastq files are automatically downloaded using the provided sample identifiers.

### Features:
(**a**) Download of fastq files either directly from ENA if not local<br/>
(**b**) Quality control of reads with [`fastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<br/>
(**c**)	Alignment of sequencing reads with [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to produce bam files.<br/>
(**d**) Filtering of bam files to remove mitochondrial reads and PCR duplicates with [`sambamba`](https://lomereiter.github.io/sambamba/index.html)<br/>
(**e**) Fragment size quality control plots with [`ATACseqQC`](https://bioconductor.org/packages/release/bioc/html/ATACseqQC.html)<br/>
(**f**) Peak calling with [`MACS2`](https://pypi.org/project/MACS2/)<br/>
(**g**) bigWig generation for IGV visualisation<br/>
(**i**) Peak annotation with [`HOMER`](http://homer.ucsd.edu/homer/)<br/>
(**j**) TF motif identification with [`Meme`](https://meme-suite.org/meme/tools/meme/)<br/>


Analyses are run in parallel and in result of error you can resume with the `-resume` parameter to re-run the pipeline starting from the previous fault.


### Testing modules
Modules can be tested using the [`pytest-workflow`](https://pypi.org/project/pytest-workflow/) framework. Module test directories within the `tests` folder contain a nextflow script and a configuration yaml file defining the test for each module.

1. Install pytest-workflow

    ```console
	conda install pytest-workflow
    ```

2. Run the tests within the project directory

    ```console
	pytest --symlink --kwdof great
    ```


