#!/bin/bash

/*
 * pipeline input parameters
 */
params.genome_file = "/loaddata2/genome/barley_chrm1.fa"
params.gtf_file = "/loaddata2/genome/barley_chrm1.gtf"
params.outdir = "./outdir"
params.reads = '/data/Comboseq-Novaseq/fastqs/*R1_001.fastq.gz'
params.miRDP2package = "/home/ubuntu/1.1.4/"
log.info """
RNASEQ-NF PIPELINE
==================
genome: ${params.genome_file}
gtf: ${params.gtf_file}
outdir: ${params.outdir}
reads:${params.reads}
miRDP2package:${params.miRDP2package}
"""
.stripIndent()
workflow {

        sampleread = Channel.fromPath(params.reads).map { file -> tuple(file.simpleName, file) }
        fastqc_ch = FASTQC(sampleread)
        MULTIQC(fastqc_ch.collect())
        trimmed = TRIM(sampleread)
        index = INDEX(params.genome_file, params.gtf_file)
        align_ch = ALIGN(index, trimmed)
        INDEXALIGNMENT(align_ch)
        trimfastqc_ch = trimFASTQC(trimmed)
        trimMULTIQC(trimfastqc_ch.collect())

        filter_ch = FILTER(align_ch, trimmed)

        GeneCount(align_ch)
        prematrix_ch= prematrix(align_ch)
        genematrix=matrix(align_ch, prematrix_ch.collect())
        mirdp_ch= miRDP2(params.genome_file, params.miRDP2package)
        mirdp_ch.view()
        miRID(trimmed, params.miRDP2package, params.genome_file, mirdp_ch.collect())
        featureCounts(params.gtf_file, trimmed, filter_ch)
        BLASTN(trimmed, index)
}

process INDEX {
cpus 8
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
cpus 8
    tag "FASTQC on ${sample_id}"
    publishDir "$params.outdir/quality.untrimmed", mode:'copy'

    input:
    tuple val(sample_id), path(sampleread)

    output:
    path "*"

    script:
    """
        mkdir ${sample_id}_fastqc
        fastqc -t 8 -o ${sample_id}_fastqc -f fastq ${sampleread}
    """
}
process MULTIQC {
cpus 8
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
cpus 8
        tag "trimming $sampleread"
        publishDir "$params.outdir/trimmed"
        input:
        tuple val(sample_id), path(sampleread)

        output:
        tuple val(sample_id), path('*')

        script:
        """
        cutadapt -u -4 -j 24 -a A{8} -o ${sample_id}.trimmed.fq.gz ${sampleread}
        """
}
process trimFASTQC {
cpus 8
    tag "FASTQC on $sample_id"
    publishDir "$params.outdir/quality.trimmed", mode:'copy'

    input:
    tuple val(sample_id), path(trimmed)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}
    fastqc -t 8 -o fastqc_${sample_id}_logs -f fastq ${trimmed}
    """
}
process trimMULTIQC {
cpus 8
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
cpus 8
    publishDir "$params.outdir/alignment", mode: 'copy', pattern:'*Aligned.sortedByCoord.out.bam'
    publishDir "$params.outdir/genecounts", mode: 'copy', pattern:'*ReadsPerGene.out.tab'
    publishDir "$params.outdir/alignmentlog", mode: 'copy', pattern:'*Log.final.out'
    input:
    path index
    tuple val(sample_id), path(trimmed)

    output:
    tuple val(sample_id), path('*')
    script:
    """
        STAR --runThreadN $task.cpus --genomeDir $index \
         --readFilesIn ${trimmed} \
         --readFilesCommand zcat \
         --sjdbInsertSave All \
         --alignIntronMax 1 \
         --outFilterMismatchNoverLmax 0.05 \
         --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
         --outFilterMatchNminOverLread 0 \
         --quantMode GeneCounts \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingThreadN  $task.cpus \
         --outFileNamePrefix $sample_id
    """
}
process INDEXALIGNMENT {
cpus 8
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
cpus 8
        publishDir "$params.outdir/separated/srna", mode: 'copy', pattern: "*srna.bam"
        publishDir "$params.outdir/separated/mrna", mode: 'copy', pattern: "*mrna.bam"
        input:
        tuple val(sample_id), path(align_ch)
        tuple val(sample_id), path(trimmed)

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
        publishDir "$params.outdir/miRDP2", mode: 'copy', pattern: "bowtie.mirdp*"
        input:
        path genome_file
        path miRDP2package

        output:
        path('bowtie.mirdp*')


        script:
        """
        echo '${params.outdir}/miRDP2/bowtie.mirdp' > indexlocation.txt
        bowtie2-build --large-index --threads 16 -f ${genome_file} bowtie.mirdp
        wget --directory-prefix=tmp.dir 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/*'
        zcat tmp.dir/* > Rfam_fa

        bowtie2-build --threads 16 -f Rfam_fa ${miRDP2package}/scripts/index/rfam_index
        """
}
process miRID {
cpus 8
publishDir "$params.outdir/miRDP2/$sample_id"

        input:
        tuple val(sample_id), path(trimmed)
        path miRDP2package
        path genome_file
        path mirdp_ch

        output:
        path "*"

        script:
        """
        fqtrim -p 8 -l 12 -C  -n read -o non-redundant.fq -o non-redundant.fq ${trimmed}
        fastq_to_fasta -Q33 -i ${sample_id}_merged.assembled.non-redundant.fq -o ${sample_id}_merged.assembled.non-redundant.fa
        bash ${miRDP2package}/miRDP2-v1.1.4_pipeline.bash -T -g ${genome_file} -x $params.outdir/miRDP2/bowtie.mirdp -f -i ${sample_id}_merged.assembled.non-redundant.fa -p 16
        cat ${mirdp_ch[0]} > test.txt
        """
}


process GeneCount {
cpus 8
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
cpus 8
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
cpus 8
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

process featureCounts {
cpus 8
        tag "Finding gene overlap"
        publishDir "$params.outdir/featurecounts/${sample_id}"
        input:
        path gtf
        tuple val(sample_id), path(trimmed)
        tuple val(type), path(filter_ch)

        output:
        path '*.txt'

        shell:
        '''
        mkdir '!{sample_id}'
        samtools view -hf 0x2 !{filter_ch[1]} |samtools view -hF 0x100 |awk  '$9>50 || $9< -50 {print $0}' |samtools sort -n |samtools view -h > '!{sample_id}.mrna.bam'
        featureCounts -p -M -B -t exon  -g gene_id -a !{gtf} -o !{sample_id}.mrna.counts.txt !{sample_id}.mrna.bam

        samtools view -hf 0x2 !{filter_ch[0]} |samtools view -hF 0x100 |awk  '$9>50 || $9< -50 {print $0}' |samtools sort -n |samtools view -h > '!{sample_id}.srna.bam'
        featureCounts -p -M -B -t exon  -g gene_id -a !{gtf} -o !{sample_id}.srna.counts.txt !{sample_id}.srna.bam
        '''
}
process BLASTN {
cpus 8
        tag "blast sequences against Lunardon et al. (2020) barley sRNA loci database to isolate sRNA"
        publishDir "$params.outdir/blastn_made_bed/${sample_id}"
        input:
        tuple val(sample_id), path(trimmed)
        path(index)

        output:
        path '*.bed'

        script:
        """
        fqtrim -l 12 -C  -n read -o non-redundant.fq  ${trimmed}
        fastq_to_fasta -Q33 -i ${sample_id}.merged.assembled.non-redundant.fq -o ${sample_id}.merged.assembled.non-redundant.fa
        makeblastdb -in ${sample_id}.merged.assembled.non-redundant.fa -dbtype nucl
        blastn -task blastn-short -db final.fasta -query /loaddata2/blastn/allloci/annotations_16768778051515.fasta -evalue 1e-6 -outfmt 6 -out ${sample_id}_blast.out
        egrep -v ^N ${sample_id}_blast.out |cut -f 2 > ${sample_id}_mappedreads.txt
        seqtk subseq ${sample_id}.merged.assembled.non-redundant.fa ${sample_id}_mappedreads.txt > ${sample_id}_out.fasta

        STAR --genomeDir ${index} -runThreadN 14 --readFilesIn ${sample_id}_out.fasta --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard --outFileNamePrefix ${sample_id}
        convert2bed -i bam -o BED < Aligned.sortedByCoord.out.bam > final.bed
        """
}
