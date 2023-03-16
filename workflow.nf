#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Optimising data
//
// Try turning on correction in fastp, corrects mismatched reads on alignment


// Variables
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
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
include {module as STAR; index as STAR_INDEX} from "./modules/star/main.nf"


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

    // Fake processes
    // Used for when Nextflow doesn't resume efficiently
    query = Channel.value("reads-trimmed")
    LOADER(
        query.combine(EXTRACT)) // Reads to trimmed reads

    reads_out = LOADER.out.reads

    emit:
    reads = reads_out
}

workflow PLANT {}
workflow FUNGI {
    take:
    reads

    main:
    id = Channel.of("${params.data.fungi}")
    genome = Channel.fromPath("${genomes}/${params.data.fungi}/genome.fasta")
    gtf = Channel.fromPath("${genomes}/${params.data.fungi}/genome.gtf")
    reads = reads

    // STAR indexer also returns genome length, under the output genome_length
    STAR_INDEX(
        id.combine(genome).combine(gtf))

    genome_length = STAR_INDEX.out.genome_length

    // STAR accepts: reads, genome, genome_length, gtf, indexer output
    STAR(
        reads
            .combine(genome)
            .combine(STAR_INDEX.out.genome_length)
            .combine(gtf)
            .combine(STAR_INDEX.out.index))

    emit:
    star = STAR.out.star
}

workflow {
    MAIN()
    FUNGI(MAIN.out.reads)
}
