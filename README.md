# GTN_sWGS
### Description
This is a pipeline to obtain relative copy number from fastq files from sWGS. This pipeline is written using nextflow. I have aligned here to hg19 for ease, as hg38 bin annotations are not easily available for analysis with QDNASeq at this time.
### Pipeline summary
#### Processes
1. `FASTQC_RAW`
	+ fastQC is used on the raw files to determine their quality
2. `TRIMMOMATIC`
	+ Trimmomatic is used to trim out adapters and polyG tails. These are available in a text file in /bin and specified with `params.adap`
3. `FASTQC_TRIM`
	+ Used to determine the quality of the sequences after trimming out adapters and polyG tails.
4. `MultiQC_TRIM`
	+ Used to collect results from fastQC before and after trim and create one report for all fastqc reports
5. `ALIGN`
	+ `bwamem` and `samtools` are used to align the sWGS files to hg19 and sort/index resulting bam files. A genome index is needed for this step and its path is specified using `params.index`. `qualimap` is also used to determine efficiency of the mapping 
6. QDNNASeq is used to obtain relative copy number from bam files using a customised Rscript in `$projectDir/bin`. This process is a little different, as two separate scripts are needed depending on whether the sex chromosomes are to be included in analysis as specified using the `--sexChrIncl` flag:
	+ `QDNASEQ` - default option. Sex chromosomes are not included in analysis, nextflow directs to the `QDNASeq.R` script in `$projectDir/bin`
	+ `QDNASEQSEX` - Sex chromosomes included in analysis, based on the presence of the `--sexChrIncl` flag in the main `nextflow run` command

#### Parameters
* `--reads` can be used to specify the path to the reads. Default is set to `$projectDir/sWGS/fastq/*_{R1,R2}.fq.gz`
* `--index` path to BWA indexed genome. Default is hg19 in reference data in ICR cluster Alma.
* `--outdir` where to publish results. Deafult `$projectDir/results`
* `--adap` path to textfile where adapter sequences are. These are essential for trimmomatic.
* `--binsize` default `30`, allows you to specify the bin size for QDNASeq. available options are 1, 5, 10, 15, 30, 50, 100, 500, and 1000 kbp.
* `--sexChrIncl` specify this flag if you want to include sex chromsomes. Othersie QDNASeq will ignore them.

### QDNASeq
Two customised scripts for QDNASeq can be found in /bin (QDNAseq.R QDNAseq_sexChrIncl.R). These have been made using the QDNASeq vignette available [here](https://bioconductor.org/packages/release/bioc/vignettes/QDNAseq/inst/doc/QDNAseq.pdf)
