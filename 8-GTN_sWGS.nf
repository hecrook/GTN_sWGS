/*
 * Gestational sWGS input parameters
 */

params.reads = "$projectDir/sWGS/fastq/*_{R1,R2}.fq.gz"
params.index = "/data/reference-data/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/"
params.outdir = "$projectDir/results"
params.adap = "$projectDir/bin/TruSeq3-PE-2.fa"
params.binsize = "30"
params.sexChrIncl = false


log.info """\
    s W G S - N F   P I P E L I N E
    ===================================
    index       : ${params.index}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    binsize      : ${params.binsize}
    """
    .stripIndent()

process FASTQC_RAW {
    conda 'fastqc'
    publishDir "$projectDir/results/fastqc/raw", mode:'copy'
    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_fastqc"

    script:
    """
    mkdir "${sample_id}_fastqc"
    fastqc -o ${sample_id}_fastqc -f fastq ${reads}
    """

}

process TRIMMOMATIC {
    publishDir "$projectDir/results/trimmomatic", mode:'copy'
    conda 'trimmomatic'
    input:
    path adap
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed_paired.fastq")

    script:
    """
    mkdir -p ${sample_id}_trimmed/

    # Trimming
    trimmomatic PE \
                ${reads[0]} ${reads[1]} \
                ${sample_id}_R1_trimmed_paired.fastq ${sample_id}_R1_trimmed_unpaired.fastq \
                ${sample_id}_R2_trimmed_paired.fastq ${sample_id}_R2_trimmed_unpaired.fastq \
                ILLUMINACLIP:$adap:2:30:10 \

"""
}

process FASTQC_TRIM {
    publishDir "$projectDir/results/fastqc/trimmed", mode:'copy'
    conda 'fastqc'
    input:
    tuple val(sample_id), path(trimreads)

    output:
    path "${sample_id}_trimmed_fastqc"

    script:
    """
    mkdir "${sample_id}_trimmed_fastqc"
    fastqc -o ${sample_id}_trimmed_fastqc -f fastq ${trimreads}
    """

}

process MULTIQC_TRIM {
    publishDir "$projectDir/results/multiqc", mode:'copy'
    conda 'multiqc'
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
    publishDir "$projectDir/results/align/hg19", mode: 'copy'
    input:
    path index
    tuple val(sample_id), path(trimreads)

    output:
    path "${sample_id}_sorted.bam", emit: bamfile
    path "${sample_id}_sorted.bam.bai"
    path "${sample_id}_qualimap_results/"

    script:
    """

    INDEX=`find -L BWAIndex/ -maxdepth 1 -name "*.amb" | sed 's/\\.amb\$//'`
    bwa mem  -M -t $task.cpus \$INDEX ${trimreads[0]} ${trimreads[1]} > ${sample_id}.sam

    samtools view -S -b ${sample_id}.sam > ${sample_id}.bam
    samtools sort ${sample_id}.bam -o ${sample_id}_sorted.bam
    samtools index ${sample_id}_sorted.bam

    mkdir -p ${sample_id}_qualimap_results

    qualimap bamqc -bam ${sample_id}_sorted.bam -outdir ${sample_id}_qualimap_results --paint-chromosome-limits --genome-gc-distr HUMAN --collect-overlap-pairs -outformat HTML

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
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
    Channel
        .fromPath(params.index, checkIfExists: true)
        .set { index_ch }
    fastqc_raw_ch = FASTQC_RAW(read_pairs_ch)
    trim_ch = TRIMMOMATIC(params.adap, read_pairs_ch) 
    fastqc_trim_ch = FASTQC_TRIM(trim_ch)
    MULTIQC_TRIM(fastqc_raw_ch.mix(fastqc_trim_ch).collect())
    align_ch = ALIGN(index_ch.first(), trim_ch)
    align_ch.bamfile.view()
    if ( params.sexChrIncl ) {
        QDNASEQSEX(align_ch.bamfile.collect(), params.binsize)
    }
    else {
        QDNASEQ(align_ch.bamfile.collect(), params.binsize)
    }
}
