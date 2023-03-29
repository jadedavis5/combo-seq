#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Optimising data
//
// Try turning on correction in fastp, corrects mismatched reads on alignment


// Variables
// reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
reads = "${params.data.path}/reads-trimmed/${params.data.glob}"
reads_unzipped = "${params.data.path}/${params.data.reads_unzipped}/${params.data.glob_unzipped}"
genomes = "${params.data.path}/${params.data.genomes}"
adapters = "${params.data.path}/${params.data.adapters}"


log.info """
# ============================= #
# - RNA-seq Nextflow Pipeline - #
# ============================= #

> reads:	${reads}
"""
.stripIndent()


// Debugging
// println(params)


// Load modules
include {module as INSPECTOR} from "./modules/inspector/main.nf"
include {module as LOADER} from "./modules/loader/main.nf"
include {module as EXTRACT} from "./modules/extract/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"
include {module as FASTQC} from "./modules/fastqc/main.nf"
include {module as MULTIQC} from "./modules/multiqc/main.nf"
include {module as TRINITY} from "./modules/trinity/main.nf"
include {length as SEQ_LENGTH} from "./modules/seqkit/main.nf"
include {align as STAR; index as STAR_INDEX} from "./modules/star/main.nf"


// Channels
// Reads
// Need to fix logic and variables here
if ("${params.workflow.skip_extract}" == "True") {
    // reads_pe = Channel.fromFilePairs(reads_unzipped)
    reads_pe = Channel.fromFilePairs(reads)
} else {
    reads_pe = Channel.fromFilePairs(reads)
}

// Adapters
adapters = Channel.fromPath(adapters)


workflow MAIN {
    main:
    // Unzip reads, unless disabled
    if ("${params.workflow.skip_extract}" == "True") {
        EXTRACT = reads_pe
    } else {
        EXTRACT(reads_pe)
        EXTRACT = EXTRACT.out.reads
    }

    // Buffer testing
    // EXTRACT = EXTRACT
    //     .combine(adapters)
    //     .buffer(size: params.fastp.buffer.toInteger(), remainder: true)
    //     .groupTuple()
    // EXTRACT.view()
    // INSPECTOR(EXTRACT).view()

    // Trimming and QC
    // FASTP(EXTRACT.combine(adapters))
    // FASTQC(FASTP.out.reads.combine(adapters))
    // MULTIQC(FASTQC.out.reports.collect())

    // LOADER
    // Used for when Nextflow doesn't resume efficiently
    // Make sure directories terminate with a trailing slash "/"
    query = Channel.value("reads-trimmed/")
    LOADER(
        query.combine(EXTRACT.flatMap {n -> n[0]})) // Reads to trimmed reads

    reads_out = LOADER.out
    // reads_out = FASTP.out.reads

    emit:
    reads = reads_out
}

workflow PLANT {
    take:
    reads

    main:
    id = Channel.of("${params.data.plant}")
    genome = Channel.fromPath("${genomes}/${params.data.plant}/genome.fasta")
    gtf = Channel.fromPath("${genomes}/${params.data.plant}/genome.gtf")
    reads = reads

    // SEQ_LENGTH(genome)
    // genome_length = SEQ_LENGTH.out.seq_length

    // STAR_INDEX(
    //     id
    //         .combine(genome)
    //         .combine(SEQ_LENGTH.out.seq_length)
    //         .combine(gtf))
    // index = STAR_INDEX.out.index

    // // STAR
    // // Returns the STAR genome index
    // query = Channel.value("${params.data.star}/")
    // LOADER(query.combine(id))
    // index = LOADER.out

    // // STAR accepts: reads, genome_length, gtf, indexer output
    // STAR(
    //     reads
    //         .combine(genome_length)
    //         .combine(gtf)
    //         .combine(index))

    // Trinity
    // Returns BAMS from STAR
    query = Channel.value("${params.data.star}/")
    LOADER(query.combine(reads.flatMap {n -> "${n[0]}/Aligned.sortedByCoord.out.bam"}))
    bam = LOADER.out
    bam.view()
    intron_max = Channel.value("${params.data.plant_intron_max}")

    // TRINITY(reads, bam, intron_max)
    // TRINITY.out[0].view()

    // emit:
    // star = STAR.out.star
}

workflow FUNGI {
    take:
    reads

    main:
    id = Channel.of("${params.data.fungi}")
    genome = Channel.fromPath("${genomes}/${params.data.fungi}/genome.fasta")
    gtf = Channel.fromPath("${genomes}/${params.data.fungi}/genome.gtf")
    reads = reads

    // STAR indexer also returns genome length, under the output genome_length
    // SEQ_LENGTH(genome)
    // STAR_INDEX(
    //     id
    //         .combine(genome)
    //         .combine(SEQ_LENGTH.out.seq_length)
    //         .combine(gtf))
    // STAR_INDEX = LOADER("wheat-cs")

    // STAR accepts: reads, genome_length, gtf, indexer output
    // Does not need the genome, only used for index generation
    // STAR(
    //     reads
    //         .combine(SEQ_LENGTH.out.seq_length)
    //         .combine(gtf)
    //         .combine(STAR_INDEX.out.index))

    // emit:
    // star = STAR.out.star
}

workflow {
    MAIN()
    // FUNGI(MAIN.out.reads)
    PLANT(MAIN.out.reads)
}
