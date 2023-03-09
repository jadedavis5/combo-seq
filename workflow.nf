#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Variables
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"
adapters = "${params.data.path}/${params.data.adapters}"


log.info """
# ============================== #
# - ComboSeq Nextflow Pipeline - #
# ============================== #

> reads:	${reads}
"""
.stripIndent()


// Debugging
// println(params)


// Load modules
// include {TEMPLATE} from "./modules/template/main.nf"
include {module as EXTRACT} from "./modules/extract/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"


// Channels
reads_pe = Channel.fromFilePairs(reads)
adapters = Channel.fromPath(adapters)


// Main outputs are under the accessor out[0], benchmarking results under time[1]
workflow {
    EXTRACT(reads_pe)
    FASTP(EXTRACT.out[0].combine(adapters))
}
