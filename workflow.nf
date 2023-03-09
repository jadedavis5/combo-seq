#!/usr/bin/env nextflow

// Nextflow configuration
nextflow.enable.dsl = 2


// Variables
reads = "${params.data.path}/${params.data.reads}/${params.data.glob}"


log.info """
# ============================== #
# - ComboSeq Nextflow Pipeline - #
# ============================== #

> reads:	${reads}
"""
.stripIndent()


// Load modules
// include {TEMPLATE} from "./modules/template/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"


// Channels
reads_pe = Channel.fromFilePairs(reads)


workflow {
    FASTP(reads_pe)
}
