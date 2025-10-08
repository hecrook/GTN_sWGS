include { samplesheetToList } from 'plugin/nf-schema'

/*
 * Gestational sWGS input parameters
 */

params.input = ""
params.genome = ""
if ( params.genome == "hg38" ) {
        params.index = "/data/reference-data/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/"
    }
else {
    if ( params.genome == "hg19" ) {
    params.index = "/data/reference-data/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/"
    }
    else {
        exit("""
        ERROR!! 
        Either no, or invalid genome has been specified. Please choose from options hg19 or hg38
        """)
    }
}
params.outdir = "$projectDir/results"
params.adap = "$projectDir/bin/TruSeq3-PE-2.fa"
params.binsize = "30"
params.sexChrIncl = false


log.info """\
    s W G S - N F   P I P E L I N E
    ===================================
    genome       : ${params.genome}
    index        : ${params.index}
    samplesheet  : ${params.input}
    outdir       : ${params.outdir}
    sexchr?      : ${params.sexChrIncl}
    binsize      : ${params.binsize}
    """
    .stripIndent()

process FASTQC_RAW {
    cpus 1
    memory '8000 MB'
    conda 'fastqc'
    publishDir "$projectDir/results/fastqc/raw", mode:'copy'
    input:
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    path "${meta.sample}_fastqc"

    script:
    """
    mkdir "${meta.sample}_fastqc"
    fastqc -o ${meta.sample}_fastqc -f fastq ${fastq1} ${fastq2}
    """

}

process TRIMMOMATIC {
    cpus 1
    memory '8000 MB'
    publishDir "$projectDir/results/trimmomatic", mode:'copy'
    conda 'trimmomatic'
    input:
    path adap
    tuple val(meta), path(fastq1), path(fastq2)

    output:
    tuple val(meta), path("${meta.sample}_R{1,2}_trimmed_paired.fastq")

    script:
    """
    mkdir -p ${meta.sample}_trimmed/

    # Trimming
    trimmomatic PE \
                ${fastq1} ${fastq2} \
                ${meta.sample}_R1_trimmed_paired.fastq ${meta.sample}_R1_trimmed_unpaired.fastq \
                ${meta.sample}_R2_trimmed_paired.fastq ${meta.sample}_R2_trimmed_unpaired.fastq \
                ILLUMINACLIP:$adap:2:30:10 \

"""
}

process FASTQC_TRIM {
    cpus 1
    memory '8000 MB'
    publishDir "$projectDir/results/fastqc/trimmed", mode:'copy'
    conda 'fastqc'
    input:
    tuple val(meta), path(trimreads)

    output:
    path "${meta.sample}_trimmed_fastqc"

    script:
    """
    mkdir "${meta.sample}_trimmed_fastqc"
    fastqc -o ${meta.sample}_trimmed_fastqc -f fastq ${trimreads}
    """

}

process MULTIQC_TRIM {
    publishDir "$projectDir/results/multiqc", mode:'copy'
    container "$projectDir/multiqc-1.20.sif"
    input:
    path "*"

    output:
    path "MultiQC_trim.html"

    script:
    """
    mkdir -p MultiQC_reports/
    multiqc . --filename MultiQC_trim
    """
}

process ALIGN {
    cpus 24
    memory '16000 MB'
    conda 'bwa samtools qualimap'
    publishDir "$projectDir/results/align/${params.genome}", mode: 'copy'
    input:
    path index
    tuple val(meta), path(trimreads)

    output:
    path "${meta.sample}_sorted.bam", emit: bamfile
    path "${meta.sample}_sorted.bam.bai"


    script:
    """

    INDEX=`find -L BWAIndex/ -maxdepth 1 -name "*.amb" | sed 's/\\.amb\$//'`
    echo -e "\nAligning ..."
    bwa mem  -M -t $task.cpus \$INDEX ${trimreads[0]} ${trimreads[1]} > ${meta.sample}.sam
    echo -e "\nConverting sam to bam ..."
    samtools view -S -b ${meta.sample}.sam > ${meta.sample}.bam
    echo -e "\nSorting bam files ..."
    samtools sort ${meta.sample}.bam -o ${meta.sample}_sorted.bam
    echo -e "\nIndexing bam files ..."
    samtools index ${meta.sample}_sorted.bam
    echo -e "\nFinished indexing"
    
    # mkdir -p ${meta.sample}_qualimap_results

    # qualimap bamqc -bam ${meta.sample}_sorted.bam -outdir ${meta.sample}_qualimap_results --paint-chromosome-limits --genome-gc-distr HUMAN --collect-overlap-pairs -outformat HTML

    """
}

process QDNASEQ {
    cpus 16
    memory '16000 MB'
    conda '/home/hcrook/.conda/envs/qdnaseq'
    publishDir "$projectDir/results/qdnaseq/${params.binsize}kb", mode: 'copy'

    input:
    path bamfiles
    val binsize

    output:
    path "*.RData"
    path "*.pdf"
    path "*.txt"


    script:
    """
    QDNAseq.R "${bamfiles}" "${binsize}"
    """
}

process QDNASEQSEX {
    cpus 16
    memory '16000 MB'
    conda '/home/hcrook/.conda/envs/qdnaseq'
    publishDir "$projectDir/results/qdnaseq/${params.binsize}kb_withSexChr", mode: 'copy'

    input:
    path bamfiles
    val binsize

    output:
    path "*.RData"
    path "*.pdf"
    path "*.txt"


    script:
    """
    QDNAseq_sexChrIncl.R "${bamfiles}" "${binsize}"
    """
}

workflow {
    Channel
        .fromList(samplesheetToList(params.input, "assets/schema_input.json"))
        .set { read_pairs_ch }
    read_pairs_ch.view()
    Channel
        .fromPath(params.index, checkIfExists: true)
        .set { index_ch }
    fastqc_raw_ch = FASTQC_RAW(read_pairs_ch)
    trim_ch = TRIMMOMATIC(params.adap, read_pairs_ch) 
    fastqc_trim_ch = FASTQC_TRIM(trim_ch)
    MULTIQC_TRIM(fastqc_raw_ch.mix(fastqc_trim_ch).collect())
    align_ch = ALIGN(index_ch.first(), trim_ch)
    align_ch.bamfile.view()
    if ( params.genome == "hg19" ) {
        if ( params.sexChrIncl ) {
            QDNASEQSEX(align_ch.bamfile.collect(), params.binsize)
        }
        else {
            QDNASEQ(align_ch.bamfile.collect(), params.binsize)
        }
    } else {}

}
