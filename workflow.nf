#!/usr/bin/env nextflow

/* Pipeline upgrade notes*/

// Nextflow configuration
nextflow.enable.dsl = 2


// Add project root to paths


log.info """
# ============================== #
# - ComboSeq Nextflow Pipeline - #
# ============================== #

> genome: ${params.genome}
> gtf: ${params.gtf}
> outdir: ${params.outdir}
> reads:${params.reads}
"""
.stripIndent()


workflow {
        read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )
        trimmed_pairs_ch = TRIM(read_pairs_ch)
        index = INDEX(params.genome, params.gtf)

        fastqc_ch = FASTQC(read_pairs_ch)
        MULTIQC(fastqc_ch.collect())
        align_ch = ALIGN(index, trimmed_pairs_ch)
        INDEXALIGNMENT(align_ch)
        trimfastqc_ch = trimFASTQC(trimmed_pairs_ch)
        trimMULTIQC(trimfastqc_ch.collect())

        filter_ch = FILTER(align_ch, trimmed_pairs_ch)

        GeneCount(align_ch)
        // prematrix_ch= prematrix(align_ch)
        // matrix(align_ch, prematrix_ch.collect())
        // mirdp_ch= miRDP2(params.genome, params.miRDP2package)
        // mirdp_ch.view()
        // miRID(trimmed_pairs_ch, params.miRDP2package, params.genome, mirdp_ch.collect())
        // featureCounts(params.gtf, trimmed_pairs_ch, filter_ch)
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
        publishDir "$params.outdir/miRDP2", mode: 'copy', pattern: "bowtie.mirdp*"
        input:
        path genome
        path miRDP2package

        output:
        path('bowtie.mirdp*')


        script:
        """
        echo '${params.outdir}/miRDP2/bowtie.mirdp' > indexlocation.txt
        bowtie-build --large-index --threads 16 -f ${genome} bowtie.mirdp
        wget --directory-prefix=tmp.dir 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF00001.fa.gz'
        zcat tmp.dir/* > Rfam_fa

        bowtie-build --threads 16 -f Rfam_fa ${miRDP2package}/scripts/index/rfam_index
        """
}

//issue with this that miRDP2 can't find mature_index, but bowtie doesn't make one (shown in script_err)
//also it says that none of them are paired
process miRID {
cpus 8
publishDir "$params.outdir/miRDP2/$sample_id"

        input:
        tuple val(sample_id), path(trimmed_pairs_ch)
        path miRDP2package
        path genome
        path mirdp_ch

        output:
        path "*"

        script:
        """
        pear -j 8 -f ${trimmed_pairs_ch[0]} -r ${trimmed_pairs_ch[1]} -j 8 -m 40 -n 12 -o ${sample_id}_merged
        fqtrim -p 8 -l 12 -C  -n read -o non-redundant.fq -o non-redundant.fq ${sample_id}_merged.assembled.fastq
        fastq_to_fasta -Q33 -i ${sample_id}_merged.assembled.non-redundant.fq -o ${sample_id}_merged.assembled.non-redundant.fa
        bash ${miRDP2package}/miRDP2-v1.1.4_pipeline.bash -T -g ${genome} -x $params.outdir/miRDP2/bowtie.mirdp -f -i ${sample_id}_merged.assembled.non-redundant.fa -p 16
        cat ${mirdp_ch[0]} > test.txt  //need to take this away or incorporate another way 
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

process featureCounts {
cpus 16
        tag "Finding gene overlap"
        publishDir "$params.outdir/featurecounts/${sample_id}"
        input:
        path gtf
        tuple val(sample_id), path(trimmed_pairs_ch)
        tuple val(type), path(filter_ch)

        output:
        path '*.txt'
        path '*.summary'

        shell:
        '''
        mkdir '!{sample_id}'
        samtools view -hf 0x2 !{filter_ch[1]} |samtools view -hF 0x100 |awk  '$9>50 || $9< -50 {print $0}' |samtools sort -n |samtools view -h > '!{sample_id}.mrna.bam'
        featureCounts -p -M -B -t exon  -g gene_id -a !{gtf} -o !{sample_id}.mrna.counts.txt !{sample_id}.mrna.bam

        samtools view -hf 0x2 !{filter_ch[0]} |samtools view -hF 0x100 |awk  '$9>50 || $9< -50 {print $0}' |samtools sort -n |samtools view -h > '!{sample_id}.srna.bam'
        featureCounts -p -M -B -t exon  -g gene_id -a !{gtf} -o !{sample_id}.srna.counts.txt !{sample_id}.srna.bam
        '''
}
