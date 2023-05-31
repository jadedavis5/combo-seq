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
include {loader1 as LOADER1; loader2 as LOADER2} from "./modules/loader/main.nf"
include {module as INSPECTOR} from "./modules/inspector/main.nf"
include {module as EXTRACT} from "./modules/extract/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"
include {module as FASTQC} from "./modules/fastqc/main.nf"
include {module as MULTIQC} from "./modules/multiqc/main.nf"
include {assemble1 as TRINITY1; assemble2 as TRINITY2} from "./modules/trinity/main.nf"
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
    path = Channel.value("reads-trimmed/")
    LOADER1(
        EXTRACT.flatMap {n -> n[0]}
            .combine(path)) // Reads to trimmed reads

    reads_out = LOADER1.out
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
    // path = Channel.value("${params.data.star}")
    // LOADER1(id.combine(path).combine(id))
    // index = LOADER1.out

    // STAR accepts: reads, genome_length, gtf, indexer output
    // STAR(
    //     reads
    //         .combine(genome_length)
    //         .combine(gtf)
    //         .combine(index))

    // // Trinity
    path = Channel.value("${params.data.star}/")
    LOADER2(reads.flatMap {n -> "${n[0]}-${params.data.plant}"}
       .combine(path)
       .combine(Channel.value("Aligned.sortedByCoord.out.bam")))
    bam = LOADER2.out
    intron_max = Channel.value("${params.data.plant_intron_max}")

    // TRINITY2: Guided assembly
    TRINITY2(reads.flatMap {n -> [["${n[0]}-${params.data.plant}"]]}
            .join(bam)
            .combine(intron_max))
    // TRINITY2.out[0].view()

    // emit:
    // star = STAR.out.star
}

workflow FUNGI {
}

workflow {
    MAIN()
    // FUNGI(MAIN.out.reads)
    PLANT(MAIN.out.reads)
}
