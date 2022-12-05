#!/bin/bash

/*
 * pipeline input parameters
 */
params.genome_file = "/data2/genome/barley_chrm1.fa"
params.gtf_file = "/data2/star_index/barley_chrm1.gtf"
params.multiqc = "/data/multiqc"
params.outdir = "/data2/outdir"
params.reads = '/data/Comboseq-Novaseq/fastqs/*R{1,2}_001.fastq.gz'
params.bam = '/data2/outdir/aligned/*{1,2}Aligned.sortedByCoord.out.bam'
log.info """
RNASEQ-NF PIPELINE
==================
genome: ${params.genome_file}
gtf: ${params.gtf_file}
outdir: ${params.outdir}
reads:${params.reads}
bam:${params.bam}
"""
.stripIndent()


workflow {
        read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
        trimmed_pairs_ch = TRIM(read_pairs_ch)
        index = INDEX(params.genome_file, params.gtf_file)
        fastqc_ch = FASTQC(read_pairs_ch)
        MULTIQC(fastqc_ch.collect())
        align_ch = ALIGN(index, trimmed_pairs_ch)
        INDEXALIGNMENT(align_ch)
        trimfastqc_ch = trimFASTQC(trimmed_pairs_ch)
        trimMULTIQC(trimfastqc_ch.collect())
        filter_ch = FILTER(align_ch, trimmed_pairs_ch)
        filter_ch.view()
        featureCounts(params.gtf_file, trimmed_pairs_ch, filter_ch)

}
process INDEX {
cpus 16
    input:
    path genome
    path gtf

    output:
    path star_index

    script:
    """
        STAR --runThreadN $task.cpus --runMode genomeGenerate --limitGenomeGenerateRAM 63000000000 --genomeDir star_index \
        --genomeFastaFiles $genome \
        --sjdbGTFfile $gtf \
        --sjdbOverhang 149 \
        --quantMode GeneCounts \
        --genomeSAindexNbases 13
   """
}
process FASTQC {
cpus 16
    tag "FASTQC on $sample_id"
    publishDir "$params.outdir/quality.untrimmed", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq ${reads}
    """
}
process MULTIQC {
cpus 16
        publishDir "$params.outdir/quality.untrimmed", mode:'copy'
        input:
        path '*'

        output:
        path 'multiqc_report.html'

        script:
        """
        multiqc .
        """
}
process TRIM {
cpus 16
        tag "trimming samples"
        publishDir "$params.outdir/trimmed"
        input:
        tuple val(sample_id), path (reads)

        output:
        tuple val(sample_id), path('*')

        script:
        """
        fqtrim -5 GUUCAGAGUUCUACAGUCCGACGAUC -3 TGGAATTCTCGGGTGCCAAGG -o trimmed.fq.gz ${reads[0]},${reads[1]}
        """
}
process trimFASTQC {
cpus 16
    tag "FASTQC on $sample_id"
    publishDir "$params.outdir/quality.trimmed", mode:'copy'

    input:
    tuple val(sample_id), path(trimmed_pairs_ch)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq ${trimmed_pairs_ch}
    """
}
process trimMULTIQC {
cpus 16
        publishDir "$params.outdir/quality.trimmed", mode:'copy'

        input:
        path '*'

        output:
        path 'multiqc_report.html'

        script:
        """
        multiqc .
        """
}

process ALIGN {
cpus 16
    publishDir "$params.outdir/alignment", mode: 'copy', pattern:'*Aligned.sortedByCoord.out.bam'


    input:
    path index
    tuple val(sample_id), path(trimmed_pairs_ch)

    output:
    path '*Aligned.sortedByCoord.out.bam'

    script:
    """
        STAR --runThreadN $task.cpus --genomeDir $index \
         --readFilesIn ${trimmed_pairs_ch[0]} ${trimmed_pairs_ch[1]} \
         --readFilesCommand zcat \
         --sjdbInsertSave All \
         --alignIntronMax 1 \
         --outFilterMismatchNoverLmax 0.05 \
         --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
         --outFilterMatchNminOverLread 0 \
         --clip3pAdapterSeq TGGAATTCTCGGGTGCCAAGG CCTTGGCACCCGAGAATTCCA --clip3pAdapterMMp 0.1 0.1 \
         --quantMode GeneCounts \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingThreadN  $task.cpus \
         --outFileNamePrefix $sample_id
    """
}
process INDEXALIGNMENT {
cpus 16
        publishDir "$params.outdir/indexed.aligment", mode: 'copy'
        input:
        path align_ch

        output:
        path "*.bai"

        script:
        """
        samtools index ${align_ch}
        """
}
process FILTER {
cpus 16
        publishDir "$params.outdir/separated/srna", mode: 'copy', pattern: "*srna.bam"
        publishDir "$params.outdir/separated/mrna", mode: 'copy', pattern: "*mrna.bam"
        input:
        path align_ch
        tuple val(sample_id), path(trimmed_pairs_ch)

        output:
        tuple val(sample_id), path('*')

        script:
        """
        samtools view -hf 2 ${align_ch} | \
        awk 'substr(\$0,1,1)=="@" ||  (\$9<=50)' | \
        samtools view -b > ${sample_id}_srna.bam

        samtools view -hf 2 ${align_ch} | \
        awk 'substr(\$0,1,1)=="@" ||  (\$9>50)' | \
        samtools view -b > ${sample_id}_mrna.bam
        """
}

process featureCounts {
cpus 16
        tag "Finding gene overlap"
        publishDir "$params.outdir/featurecounts"
        input:
        path gtf
        tuple val(sample_id), path(trimmed_pairs_ch)
        tuple val(type), path(filter_ch)

        output:
        path '*.txt'
        path '*.summary'

        script:
        """
        featureCounts -T 4 -s 2 -a $gtf -o ${sample_id}.mrna.featurecounts.txt ${filter_ch[1]}
        featureCounts -T 4 -s 2 -a $gtf -o ${sample_id}.srna.featurecounts.txt ${filter_ch[0]}
        """
}
