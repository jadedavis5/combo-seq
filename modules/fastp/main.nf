// Module information
name = "fastp"
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
if (module.memory == "-1") {
    memory = Math.floor(module.cores.toInteger() * 0.5)
    module.memory = "${memory}G"
    binding.variables.remove 'memory'
} else if (module.memory == "0") {
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
flags = ""
if (module.autodetect == "True") {
    flags += " --detect_adapter_for_pe"
}
flags = flags.strip()

chunks = (module.cores).toInteger() * 2


process module {
    tag "${name}-${id}"
    cpus module.cores
    memory module.memory
    time module.time
    queue module.queue

    publishDir "${params.data.out}/${params.data.reads_trimmed}",
        mode: "copy",
        overwrite: true,
        pattern: "trim-${id}*"
    publishDir "${params.data.reports}/${params.data.fastp}",
        mode: "copy",
        overwrite: true,
        pattern: "fastp*",
        saveAs: {"${id}-${it}"}
    publishDir "${params.data.out}/${params.data.time}",
        mode: "copy",
        overwrite: true,
        pattern: "time-${id}-${name}.txt"

    input:
    tuple val(id), path(reads), path(adapters)

    output:
    tuple val(id), path("trim-${id}*"), emit: "reads"
    path("fastp*"), emit: "reports"
    path("time-${id}-${name}.txt"), emit: "tfile", optional: true

    shell:
    '''
    flags="!{module.flags}"
    if [[ "!{module.autodetect}" == "True" ]]; then
        flags+="${flags} !{flags}"
    else
        flags+="${flags} --adapter_fasta !{adapters[0]}"
    fi

    if [[ "!{params.workflow.compress}" == "True" ]] && [[ "!{module.compress}" == "True" ]]; then
        flags+="${flags} -z !{module.compress_level}"
    fi

    tfile="time-!{id}-!{name}.txt"
    !{params.time.bin} !{params.time.flags} -o "${tfile}" \
        !{bin} ${flags} \
        -i "!{reads[0]}" -I "!{reads[1]}" \
        -f "!{module.trimf}" -F "!{module.trimf}" \
        -t "!{module.trimt}" -T "!{module.trimt}" \
        -o "trim-!{reads[0]}" -O "trim-!{reads[1]}" \
        -q "!{module.quality}"

    # Process information
    echo "	Allocated resources" >> "${tfile}"
    echo "	CPUs: !{task.cpus}" >> "${tfile}"
    echo "	Threads: !{threads}" >> "${tfile}"
    echo "	Memory: !{task.memory}" >> "${tfile}"
    echo "	Time: !{task.time}" >> "${tfile}"
    echo "	Queue: !{task.queue}" >> "${tfile}"
    '''
}
