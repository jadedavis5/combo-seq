#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


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
include {module as EXTRACT} from "./modules/extract/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"


// Channels
reads_pe = Channel.fromFilePairs(reads)
adapters = Channel.fromPath(adapters)


// Main outputs are under the accessor out[0], benchmarking results under time[1]
workflow {
    // Unzip reads, unless disabled
    if ("${params.workflow.skip_extract}" == "True") {
        EXTRACT = reads_pe
    } else {
        EXTRACT(reads_pe)
        EXTRACT = EXTRACT.out[0]
    }

    // Trimming and QC
    FASTP(EXTRACT.combine(adapters))

    // Index genome using STAR
}

workflow plant {}
workflow fungi {}
