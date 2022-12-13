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
params.miRDP2 = "/home/ubuntu/1.1.4/miRDP2-v1.1.4_pipeline.bash"
log.info """
RNASEQ-NF PIPELINE
==================
genome: ${params.genome_file}
gtf: ${params.gtf_file}
outdir: ${params.outdir}
reads:${params.reads}
bam:${params.bam}
miRDP2:${params.miRDP2}
"""
.stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */

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

        GeneCount(align_ch)
        prematrix_ch= prematrix(align_ch)
        matrix(align_ch, prematrix_ch.collect())
        miRDP2(trimmed_pairs_ch, params.miRDP2, params.genome_file)

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
        tag "trimming $sample_id"
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
    publishDir "$params.outdir/genecounts", mode: 'copy', pattern:'*ReadsPerGene.out.tab'

    input:
    path index
    tuple val(sample_id), path(trimmed_pairs_ch)

    output:
    tuple val(sample_id), path('*')
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
        tuple val(sample_id), path(align_ch)

        output:
        path "*.bai"

        script:
        """
        samtools index ${align_ch[0]}
        """
}

process FILTER {
cpus 16
        publishDir "$params.outdir/separated/srna", mode: 'copy', pattern: "*srna.bam"
        publishDir "$params.outdir/separated/mrna", mode: 'copy', pattern: "*mrna.bam"
        input:
        tuple val(sample_id), path(align_ch)
        tuple val(sample_id), path(trimmed_pairs_ch)

        output:
        tuple val(sample_id), path('*')

        script:
        """
        samtools view -hf 2 ${align_ch[0]} | \
        awk 'substr(\$0,1,1)=="@" ||  (\$9<=50)' | \
        samtools view -b > ${sample_id}_srna.bam

        samtools view -hf 2 ${align_ch[0]} | \
        awk 'substr(\$0,1,1)=="@" ||  (\$9>50)' | \
        samtools view -b > ${sample_id}_mrna.bam
        """
}

process miRDP2 {
cpus 16
        publishDir "$params.outdir/miRDP2/$sample_id"
        input:
        tuple val(sample_id), path(trimmed_pairs_ch)
        path miRDP2
        path genome_file

        output:
        path "*"

        script:
        """
        bowtie-build --large-index --threads 16 -f ${genome_file} barley.genome
        wget --directory-prefix=tmp.dir 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/*'
        zcat tmp.dir/* > Rfam.fa

        bowtie-build --threads 16 -f Rfam.fa rfam_index
        mv rfam_index.* ~/1.1.4/scripts/index/

        pear -f ${trimmed_pairs_ch[0]} -r ${trimmed_pairs_ch[1]} -j 8 -m 40 -n 12 -o ${sample_id}_merged
        fqtrim -l 12 -C  -n read -o non-redundant.fq -o non-redundant.fq ${sample_id}_merged.assembled.fastq
        fastq_to_fasta -Q33 -i ${sample_id}_merged.assembled.non-redundant.fq -o merged.assembled.non-redundant.fa
        mkdir -p "$params.outdir/miRDP2/${sample_id}/_merged.assembled.non-redundant"
        bash ${miRDP2} -g ${genome_file} -x barley.genome-f -i ${sample_id}_merged.assembled.non-redundant.fa -p 16
        """
}
process GeneCount {
cpus 16
        publishDir "$params.outdir/mapped.overall"
        input:
        tuple val(sample_id), path(align_ch)

        output:
        path '*'
        script:
        """
        echo "total reads overlapping" > '${sample_id}_overall.reads.mapped.txt'
        perl -lane '\$s+=\$F[3] ;END {print \$s}' ${align_ch[4]} >> '${sample_id}_overall.reads.mapped.txt'
        """
}
process prematrix{
cpus 16
        publishDir "$params.outdir/genematrix/tmp", mode: 'copy', pattern: "*_field4.txt"
        input:
        tuple val(sample_id), path(align_ch)

        output:
        path('*_field4.txt')

        shell:
        '''
        echo !{sample_id} > '!{sample_id}_field4.txt'
        egrep -v ^N !{align_ch[4]} |cut -f 4 >> !{sample_id}_field4.txt
        '''
}


process matrix {
cpus 16
        publishDir "$params.outdir/genematrix"
        input:
        tuple val(sample_id), path(align_ch)
        file(prematrix_ch)

        output:
        path('genematrix.txt')

        shell:
        '''
        echo gene > 'genes.txt'
        egrep -v ^N !{align_ch[4]} |cut -f 1 >> 'genes.txt'

        ls !{prematrix_ch} > 'field4.txt'
        paste 'genes.txt' $(printf "%s " $(cat 'field4.txt')) > 'genematrix.txt'
        '''
}
