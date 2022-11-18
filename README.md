
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

### Analyse an example dataset

Try the pipeline on an example dataset (all inputs will be automatically downloaded): -

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)

3. [`Configure`](https://www.nextflow.io/docs/latest/config.html) the resource profile for your HPC or local computer. A template for slurm schedulers is provided as an example in `nextflow.config`

4. Download the pipeline and test on the example dataset:

    ```console
     nextflow clone CBFLivUni/SkeletalVis-ATACSeq
     cd SkeletalVis-ATACSeq
     nextflow run main.nf -profile slurm -params-file params/example.yaml
     ```
### Analyse your own data

1. Define the sampleTable

Create a tab seperated table with unique Sample names, SRR accession numbers (if download is needed) and any additional metadata e.g

|Sample|File|Condition|
| ---|---|---|
|Control_1|SRRXXX	|Control|
|Control_2|	SRRXXX	|Control|
|Treated_1|	SRRXXX	|Treated|
|Treated_2|	SRRXXX	|Treated|

2. Define the configuration

Most parameters are set to sensible defaults within the main nextflow script, with only 5 parameters required to be altered with typical use:

|Parameter|Description|Options|
| ---|---|---|
|accession|The GEO accession of the data - used to name output data and download fastq files||
|species|The species the reads originate from - used to create the bowtie2 index	|human, mouse|
|comparison|The column in the sample table that defines the contrast to make with diffBind|e.g Treatment, Condition|
|Numerator|The column in the sample table that defines the baseline in the contrast to make with diffBind|e.g Control, Wildtype|
|Denominator|The column in the sample table that defines the denominator in the comparison with diffBind|e.g Disease, Treated|

Parameters should be defined within a yaml file. See `params/example.yaml` for an example.

3. Run the pipeline with your own parameters

    ```console
     nextflow run main.nf -profile slurm -params-file ownData.yaml
    ```

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


