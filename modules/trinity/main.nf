// Module information
name = "trinity"
module = params[name]


// Enable containers and set binary
if (module.container == "True") {
    // Check if monolithic containers are enabled, and set container
    if (params.container.monolithic == "True") {
        container = params.container.path
    } else {
        container = module.path
    }

    // Set executable
    if (module.direct == "True") {
        bin = "${container} ${module.bin}"
    } else {
        bin = "${params.container.bin} ${params.container.opts} \
        ${params.container.exec} ${params.container.exec_opts} \
        ${container} ${module.bin}"
    }
} else {
    bin = "${module.bin}"
}


// Set threads
if (params.hardware.smt == "True") {
    threads = (module.cores.toInteger() * 2).toInteger()
}


// Set memory
if (module.memory == -1) {
    throw new Exception(
        "Module ${name} does not support dynamic memory assignment, \
        please set to 0 or a positive integer")
} else if (module.memory == 0) {
    memory = 1
    module.memory = "${module.memory}G"
    binding.variables.remove 'memory'
}
else {
    module.memory = "${module.memory}G"
}


// Set time
if (!module.time) {
    println("Setting non-existent time value to 0")
    module.time_align = 0
    module.time_index = 0
}


// Module settings
logfile = ""


process assemble1 {
    tag "${name}-${id}"
    cpus module.cores
    memory module.memory
    time "${module.time}.hour"
    queue module.queue

    publishDir "${params.data.out}/${params.data.trinity}",
        mode: "copy",
        overwrite: true,
        pattern: "trinity_out_dir.Trinity.fasta",
        saveAs: {"${id}.fasta"}
    publishDir "${params.data.out}/${params.data.trinity}",
        mode: "copy",
        overwrite: true,
        pattern: "trinity_out_dir.Trinity.fasta.gene_trans_map",
        saveAs: {"${id}.gene_trans_map"}
    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"
    publishDir "${params.data.logs}/${name}",
        mode: "copy",
        overwrite: true,
        pattern: "${logfile}"

    input:
    tuple val(id), path(reads), val(intron_max)

    output:
    stdout
    tuple val(id), path("trinity_out_dir.Trinity.fasta{,.gene_trans_map}"), emit: "${name}"
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true
    path("${logfile}"), emit: "logs", optional: true

    shell:
    '''
    # flags="!{module.flags}"
    flags=""
    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} \
        --seqType fq --left "!{reads[0]}" --right "!{reads[1]}" \
        --CPU !{module.cores} \
        --max_memory !{module.memory}

    # Remove temporary work folder
    rm -rf trinity_out_dir

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}

process assemble2 {
    tag "${name}-${id}"
    cpus module.cores
    memory module.memory
    time "${module.time}.hour"
    queue module.queue

    publishDir "${params.data.out}/${params.data.trinity}",
        mode: "copy",
        overwrite: true,
        pattern: "${id}.fasta"
    publishDir "${params.data.out}/${params.data.trinity}",
        mode: "copy",
        overwrite: true,
        pattern: "${id}.gene_trans_map"
    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"
    publishDir "${params.data.logs}/${name}",
        mode: "copy",
        overwrite: true,
        pattern: "${logfile}"

    input:
    tuple val(id), path(bam), val(intron_max)

    output:
    stdout
    tuple val(id), path("${id}.{fasta,gene_trans_map}"), emit: "${name}"
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true
    path("${logfile}"), emit: "logs", optional: true

    shell:
    '''
    # flags="!{module.flags}"
    flags=""
    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} \
        --genome_guided_bam "!{bam}" \
        --genome_guided_max_intron !{intron_max} \
        --CPU !{threads} \
        --max_memory !{module.memory}

    # Move Trinity outputs
    mv trinity_out_dir/Trinity-GG.fasta "!{id}.fasta"
    mv trinity_out_dir/Trinity-GG.fasta.gene_trans_map "!{id}.gene_trans_map"

    # Remove temporary work folder
    rm -rf trinity_out_dir

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Threads: !{threads}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
