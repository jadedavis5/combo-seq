#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Optimising data
//
// Try turning on correction in fastp, corrects mismatched reads on alignment


// Variables
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
reads_unzipped = "${params.data.path}/${params.data.reads_unzipped}/${params.data.glob_unzipped}"
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
include {module as EXTRACT} from "./modules/extract/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"
include {module as FASTQC} from "./modules/fastqc/main.nf"
include {module as MULTIQC} from "./modules/multiqc/main.nf"


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


workflow {
    // Unzip reads, unless disabled
    if ("${params.workflow.skip_extract}" == "True") {
        EXTRACT = reads_pe
    } else {
        EXTRACT(reads_pe)
        EXTRACT = EXTRACT.out[0]
    }

    // Buffer testing
    // EXTRACT = EXTRACT
    //     .combine(adapters)
    //     .buffer(size: params.fastp.buffer.toInteger(), remainder: true)
    //     .groupTuple()
    // EXTRACT.vi
    // INSPECTOR(EXTRACT).view()

    // Trimming and QC
    FASTP(EXTRACT.combine(adapters))
    // FASTQC(FASTP.out.reads.combine(adapters))
    // MULTIQC(FASTQC.out.reports.collect())

    // Index genome using STAR
}

workflow plant {}
workflow fungi {}
