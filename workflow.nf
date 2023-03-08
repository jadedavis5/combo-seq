#!/usr/bin/env nextflow

/* Pipeline upgrade notes*/

// Nextflow configuration
nextflow.enable.dsl = 2


log.info """
# ============================== #
# - ComboSeq Nextflow Pipeline - #
# ============================== #

> genome:	${params.data.genome}
> gtf:  	${params.data.gtf}
> out:  	${params.data.out}
> reads:	${params.data.reads}/${params.data.glob}
"""
.stripIndent()


// Load modules
// include {TEMPLATE} from "./modules/template/main.nf"
include {module as FASTP} from "./modules/fastp/main.nf"


// Channels
reads_pe = Channel.fromFilePairs("${params.data.reads}/${params.data.glob}")


workflow {
    FASTP(reads_pe)
}
